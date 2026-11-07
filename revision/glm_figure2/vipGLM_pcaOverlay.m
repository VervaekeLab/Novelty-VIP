function vipGLM_pcaOverlay(pathName, files, cellTableIn)
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_PCAOVERLAY  Reproduce the manuscript PCA scatter of the k-means
% clustering, then overlay GLM-classified cells.
%
% INPUT  interactive: STEP 1 COHORT RESULTS file (ONE .mat); STEP 2 SESSION DATA files (one or more).
% Call with no args -> interactive (prompts for results + session files).
% Call vipGLM_pcaOverlay(pathName, files, cellTable) -> headless (used by the
% cohort script to produce the overlay at the end of a run).

    rng(1);
    headless = nargin>=3 && ~isempty(cellTableIn);
    if headless
        T=cellTableIn; T=T(T.context=="Familiar" & T.ok,:);
        sp=pathName; sf=files; if ischar(sf), sf={sf}; end
    else
        % ---- STEP 1: GLM results ----
        uiwait(msgbox({'STEP 1 of 2:','','Select the SAVED COHORT RESULTS file', ...
            '(cohortGLM_*.mat with res.cellTable).'},'Step 1','help','modal'));
            uiwait(msgbox({'Select the COHORT RESULTS file.','','ONE file: the cohortGLM_*.mat saved by vipGLM_cohort (has res.cellTable).','','NOT a session/m8xxx file.'},'Cohort results (ONE file)','help'));
        [rf,rp]=uigetfile('*.mat','COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes);
        if isequal(rf,0), return; end
        R=load(fullfile(rp,rf));
        if ~isfield(R,'res')||~isfield(R.res,'cellTable'), error('No res.cellTable.'); end
        T=R.res.cellTable; T=T(T.context=="Familiar" & T.ok,:);
        % ---- STEP 2: sessions ----
        uiwait(msgbox({'STEP 2 of 2:','','Select the SESSION .mat files (multi-select).', ...
            'These provide binnedRoisDff for the PCA.'},'Step 2','help','modal'));
            uiwait(msgbox({'Select the SESSION DATA files.','','ONE OR MORE m8xxx session .mat files (hold Cmd/Shift for several).'},'Session data (one or more)','help'));
        [sf,sp]=uigetfile('*.mat','SESSION DATA files (ONE OR MORE m8xxx .mat; multi-select)',dfSess,'MultiSelect','on');
        if isequal(sf,0), return; end
        if ischar(sf), sf={sf}; end
    end
    hasVelSign = any(strcmpi(T.Properties.VariableNames,'velSign'));

    % ---- build clustering input EXACTLY as vipGLM_clusterVsGLM / lab code ----
    binnedFiringDff=[]; cellSession=strings(0); cellIDvec=[]; velSignAll=[];
    for si=1:numel(sf)
        fn=sf{si}; S=load(fullfile(sp,fn));
        if ~isfield(S,'sData'), warning('skip %s',fn); continue; end
        sData=S.sData;
        if ~isfield(sData.imdata,'binnedRoisDff'), warning('skip %s (no binnedRoisDff)',fn); continue; end
        B=sData.imdata.binnedRoisDff; sz=size(B);
        Bre=reshape(permute(B,[2 1 3]),[sz(1)*sz(2),sz(3)]);
        binnedFiringDff=[binnedFiringDff, Bre]; %#ok<AGROW>
        nC=sz(3);
        cellSession=[cellSession; repmat(string(fn),nC,1)]; %#ok<AGROW>
        cellIDvec=[cellIDvec; (1:nC)']; %#ok<AGROW>
        % velocity sign fallback
        vsVec=nan(nC,1);
        if isfield(sData.imdata,'roiSignals')&&numel(sData.imdata.roiSignals)>=2 ...
                && isfield(sData.behavior.signals,'velocityDs') ...
                && isfield(sData.behavior.signals,'trialIndexesDs')
            dff=sData.imdata.roiSignals(2).dff; nF=size(dff,2);
            vel=double(sData.behavior.signals.velocityDs(:));
            ti =double(sData.behavior.signals.trialIndexesDs(:));
            Nf=min([nF,numel(vel),numel(ti)]); vel=vel(1:Nf); ti=ti(1:Nf);
            ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2);
            ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
            vtr=ti>0&ti<=length(ctx); cc=nan(Nf,1); cc(vtr)=ctx(ti(vtr)); cc=fillmissing(cc,'nearest');
            fam=(cc==ctx(1));
            for c=1:min(nC,size(dff,1))
                y=double(dff(c,1:Nf))'; vsVec(c)=signLocal(vel(fam),y(fam));
            end
        end
        velSignAll=[velSignAll; vsVec]; %#ok<AGROW>
    end
    if isempty(binnedFiringDff), error('No binnedRoisDff found.'); end

    data=binnedFiringDff(:,:)';
    event_time=125*30+1; binIndex=1:size(binnedFiringDff,1);
    B1_idx=find(binIndex>=30 & binIndex<event_time);
    wd=data(:,B1_idx);
    wz=(wd-mean(wd,2))./std(wd,0,2); wz(isnan(wz))=0;
    [~,score]=pca(wz);
    nPCs=min(3,size(score,2));
    clust=kmeans(score(:,1:nPCs),3,'Replicates',10);

    PC=score(:,1:3); nCellTot=size(PC,1);

    % cluster -> color matching the manuscript pie: biggest=yellow (40%),
    % middle=red (32%), smallest=blue (28%).
    red=[0.85 0.33 0.10]; yellow=[0.93 0.69 0.13]; blue=[0.00 0.45 0.74];
    bySize=[yellow; red; blue];              % order: biggest, middle, smallest
    sizes=arrayfun(@(k)sum(clust==k),1:3);
    [~,ord]=sort(sizes,'descend');           % ord(1)=biggest cluster
    colorOfClust=zeros(3,3); for r=1:3, colorOfClust(ord(r),:)=bySize(r,:); end
    fprintf('Cluster sizes: %s  (biggest->yellow, middle->red, smallest->blue)\n', mat2str(sizes));

    % ---- map GLM labels onto the clustered cells ----
    glmLabel=strings(nCellTot,1); glmVel=nan(nCellTot,1);
    for i=1:nCellTot
        rr=T(strcmpi(string(T.session),cellSession(i)) & T.cellID==cellIDvec(i),:);
        if isempty(rr), glmLabel(i)="none"; continue; end
        glmLabel(i)=string(rr.primaryFeature(1));
        if hasVelSign && ~isnan(rr.velSign(1)), glmVel(i)=sign(rr.velSign(1));
        else, glmVel(i)=velSignAll(i); end
    end
    isPupil = glmLabel=="Pupil";
    isVelP  = glmLabel=="Velocity" & glmVel>=0;
    isVelM  = glmLabel=="Velocity" & glmVel<0;

    fprintf('Clustered %d cells. GLM Pupil=%d, Vel+=%d, Vel-=%d\n', ...
        nCellTot, sum(isPupil), sum(isVelP), sum(isVelM));

    % ---- figure: 3 panels ----
    overlays={isPupil,isVelP,isVelM};
    titles={'GLM Pupil','GLM Velocity (+)','GLM Velocity (-)'};
    hiFace=[0.05 0.80 0.85; 0.10 0.95 0.20; 0.90 0.10 0.85];  % cyan, green, magenta
    figure('Color','w','Position',[40 80 1500 480]);
    for p=1:3
        subplot(1,3,p); hold on;
        % per-cluster convex-hull outline (PC1-PC2 top view) + faded points
        for k=1:3
            sel=(clust==k); P2=PC(sel,1:2);
            % hull outline
            if sum(sel)>=3
                try
                    kk=convhull(P2(:,1),P2(:,2));
                    plot(P2(kk,1),P2(kk,2),'-','Color',colorOfClust(k,:),'LineWidth',2);
                    fill(P2(kk,1),P2(kk,2),colorOfClust(k,:),'FaceAlpha',0.06,'EdgeColor','none');
                catch, end
            end
            scatter(P2(:,1),P2(:,2),26,colorOfClust(k,:),'filled', ...
                'MarkerFaceAlpha',0.30,'MarkerEdgeColor','none');
        end
        % overlay GLM group: large, black edge + bright face
        sel=overlays{p};
        scatter(PC(sel,1),PC(sel,2),80,hiFace(p,:),'filled', ...
            'MarkerEdgeColor','k','LineWidth',0.75);
        xlabel('PC1'); ylabel('PC2');
        title(sprintf('%s (n=%d) on cluster PCA',titles{p},sum(sel)));
        grid on; axis tight; box off; set(gca,'TickDir','out');
    end
    sgtitle({'GLM groups overlaid on k-means cluster PCA (familiar B1, PC1-PC2 top view)', ...
        'background + outline = cluster (red/yellow/blue); bright = GLM group'}, ...
        'FontWeight','bold');
end

function s = signLocal(x, y)
    n=min(length(x),length(y)); r=corr(x(1:n),y(1:n),'rows','complete');
    if isnan(r), s=1; else, s=sign(r); if s==0, s=1; end, end
end
