function vipGLM_clusterVsGLM()
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_CLUSTERVSGLM  Do the k-means functional clusters (Fig.1 i/j) align
% with the GLM feature groups (Vel+, Vel-, Pupil)?
%
% INPUT  STEP 1: COHORT RESULTS file (ONE .mat); STEP 2: SESSION DATA files (one or more).
% OUTPUT prints/plots cluster-vs-GLM correspondence.
%
%   Reproduces the manuscript clustering: PCA on position-binned activity of
%   ALL trials in the FAMILIAR environment (cells pooled across sessions),
%   then k-means (k=3) on the first 3 principal components. Then pulls the GLM
%   winning-feature labels from a saved cohort cellTable, matches cells by
%   session+cellID, and draws three Venn diagrams:
%       pupil/arousal cluster  vs  GLM-Pupil cells
%       velocity cluster A     vs  GLM-Vel(+) cells
%       velocity cluster B     vs  GLM-Vel(-) cells
%
%   The pupil/arousal cluster is identified by its session-ONSET TRANSIENT
%   (high activity in the first ~10 familiar trials) -- the manuscript's
%   defining property of Cluster 1 -- NOT by k-means' arbitrary cluster number.
%
%   Prompts for (1) the saved cohort results .mat (res.cellTable) and
%   (2) the session .mat files (for activity + velocity sign).

    rng(1);   % reproducible k-means

    % ---- inputs ----
    uiwait(msgbox({'STEP 1 of 2:','','Select the SAVED COHORT RESULTS file', ...
        '(cohortGLM_*.mat with res.cellTable).'},'Step 1','help','modal'));
        uiwait(msgbox({'Select the COHORT RESULTS file.','','ONE file: the cohortGLM_*.mat saved by vipGLM_cohort (has res.cellTable).','','NOT a session/m8xxx file.'},'Cohort results (ONE file)','help'));
    [rf,rp]=uigetfile('*.mat','COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes);
    if isequal(rf,0), return; end
    R=load(fullfile(rp,rf));
    if ~isfield(R,'res')||~isfield(R.res,'cellTable'), error('No res.cellTable.'); end
    T=R.res.cellTable; T=T(T.context=="Familiar" & T.ok,:);

    uiwait(msgbox({'STEP 2 of 2:','','Select the SESSION .mat files (multi-select).', ...
        'These provide the familiar-block activity for clustering.'},'Step 2','help','modal'));
        uiwait(msgbox({'Select the SESSION DATA files.','','ONE OR MORE m8xxx session .mat files (hold Cmd/Shift for several).'},'Session data (one or more)','help'));
    [sf,sp]=uigetfile('*.mat','SESSION DATA files (ONE OR MORE m8xxx .mat; multi-select)',dfSess,'MultiSelect','on');
    if isequal(sf,0), return; end
    if ischar(sf), sf={sf}; end

    % ---- build the clustering input EXACTLY as the lab's Fig.1 code does ----
    % Per session: reshape binnedRoisDff via permute([2 1 3]) then concatenate
    % horizontally across sessions. Track (session,cellID) per cell so we can
    % match cluster labels to GLM labels afterwards.
    binnedFiringDff = [];
    cellSession = strings(0);   % session filename per clustered cell (in order)
    cellIDvec   = [];           % cellID per clustered cell (column order within session)
    velSignAll  = [];           % per-cell velocity sign recomputed from session (fallback)
    for si=1:numel(sf)
        fn=sf{si}; S=load(fullfile(sp,fn));
        if ~isfield(S,'sData'), warning('skip %s',fn); continue; end
        sData=S.sData;
        if ~isfield(sData.imdata,'binnedRoisDff')
            warning('skip %s (no imdata.binnedRoisDff)',fn); continue;
        end
        B = sData.imdata.binnedRoisDff;
        sz = size(B);                                   % [d1 d2 d3] as stored
        Bre = reshape(permute(B,[2 1 3]), [sz(1)*sz(2), sz(3)]);  % literal student reshape
        binnedFiringDff = [binnedFiringDff, Bre]; %#ok<AGROW>
        nC = sz(3);
        cellSession = [cellSession; repmat(string(fn),nC,1)]; %#ok<AGROW>
        cellIDvec   = [cellIDvec; (1:nC)']; %#ok<AGROW>
        % recompute velocity sign per cell (familiar block), as a fallback for
        % results files saved before velSign was added to the cellTable
        vsVec = nan(nC,1);
        if isfield(sData.imdata,'roiSignals') && numel(sData.imdata.roiSignals)>=2 ...
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
                y=double(dff(c,1:Nf))';
                vsVec(c)=signLocal(vel(fam), y(fam));
            end
        end
        velSignAll = [velSignAll; vsVec]; %#ok<AGROW>
    end
    if isempty(binnedFiringDff), error('No binnedRoisDff found in selected sessions.'); end

    data = binnedFiringDff(:,:)';                        % [cells x (bins*trials)]
    % windows exactly as in the lab code
    event_time = 125*30+1; binIndex = 1:size(binnedFiringDff,1);
    B1_idx = find(binIndex >= 30 & binIndex < event_time);   % familiar block, drop first 29 bins
    window_idx = B1_idx;
    window_data = data(:, window_idx);
    % z-score per cell across the window
    window_data_z = (window_data - mean(window_data,2)) ./ std(window_data,0,2);
    window_data_z(isnan(window_data_z)) = 0;

    % PCA -> top 3 PCs -> k-means (k=3, 10 replicates), as in the lab code
    [~, score] = pca(window_data_z);
    nPCs = min(3,size(score,2));
    [clust, ~] = kmeans(score(:,1:nPCs), 3, 'Replicates', 10);

    % ---- onset-transient score per clustered cell (to identify the arousal cluster) ----
    % first 10 familiar trials (bins 1..125*10) minus the rest of B1, per cell
    onsetBins = find(binIndex>=30 & binIndex<125*10+1);
    restBins  = find(binIndex>=125*10+1 & binIndex<event_time);
    onset = mean(data(:,onsetBins),2,'omitnan') - mean(data(:,restBins),2,'omitnan');

    % ---- match each clustered cell to its GLM label (by session + cellID) ----
    nCellTot = size(data,1);
    glmLabel = strings(nCellTot,1); glmVelSign = nan(nCellTot,1);
    for i=1:nCellTot
        rr = T(strcmpi(string(T.session),cellSession(i)) & T.cellID==cellIDvec(i), :);
        if isempty(rr), glmLabel(i)="none"; continue; end
        glmLabel(i)=string(rr.primaryFeature(1));
        % prefer stored velocity sign; else use the value recomputed from session
        if any(strcmpi(rr.Properties.VariableNames,'velSign')) && ~isnan(rr.velSign(1))
            glmVelSign(i)=sign(rr.velSign(1));
        else
            glmVelSign(i)=velSignAll(i);
        end
    end
    fprintf('Clustered %d cells across %d sessions.\n', nCellTot, numel(sf));

    % ---- PCA + k-means done above; proceed to identify clusters ----
    onsetByClust=arrayfun(@(k) mean(onset(clust==k),'omitnan'), 1:3);
    [~,pupilClust]=max(onsetByClust);
    % the other two clusters -> match to Vel(+)/Vel(-) by velocity-sign composition
    others=setdiff(1:3,pupilClust);
    % mean velocity sign of GLM-velocity cells within each other cluster
    velSignByClust=nan(1,3);
    for k=others
        vsk=glmVelSign(clust==k & glmLabel=="Velocity");
        velSignByClust(k)=mean(vsk,'omitnan');
    end
    % assign: cluster with more positive mean sign -> Vel(+), other -> Vel(-)
    if isnan(velSignByClust(others(1))) && isnan(velSignByClust(others(2)))
        velPlusClust=others(1); velMinusClust=others(2);   % fallback
    elseif nan2neg(velSignByClust(others(1))) >= nan2neg(velSignByClust(others(2)))
        velPlusClust=others(1); velMinusClust=others(2);
    else
        velPlusClust=others(2); velMinusClust=others(1);
    end

    % ---- define the cell SETS for the Venns ----
    isClustPupil = (clust==pupilClust);
    isClustVelP  = (clust==velPlusClust);
    isClustVelM  = (clust==velMinusClust);
    isGLMpupil   = (glmLabel=="Pupil");
    isGLMvelP    = (glmLabel=="Velocity" & glmVelSign>0);
    isGLMvelM    = (glmLabel=="Velocity" & glmVelSign<0);

    % ---- three Venn diagrams ----
    figure('Color','w','Position',[80 150 1200 420]);
    subplot(1,3,1); drawVenn(isClustPupil, isGLMpupil, ...
        sprintf('Arousal cluster (n=%d)',sum(isClustPupil)), ...
        sprintf('GLM Pupil (n=%d)',sum(isGLMpupil)), [0.64 0.08 0.18]);
    title('Pupil / arousal');
    subplot(1,3,2); drawVenn(isClustVelP, isGLMvelP, ...
        sprintf('Vel cluster A (n=%d)',sum(isClustVelP)), ...
        sprintf('GLM Vel(+) (n=%d)',sum(isGLMvelP)), [0.85 0.33 0.10]);
    title('Velocity (+)');
    subplot(1,3,3); drawVenn(isClustVelM, isGLMvelM, ...
        sprintf('Vel cluster B (n=%d)',sum(isClustVelM)), ...
        sprintf('GLM Vel(-) (n=%d)',sum(isGLMvelM)), [0.93 0.69 0.13]);
    title('Velocity (-)');
    sgtitle('k-means cluster vs GLM group overlap (familiar block)','FontWeight','bold');

    % ---- console summary ----
    fprintf('\nCluster sizes: %s  (pupil/arousal cluster = #%d)\n', mat2str(histcounts(clust,1:4)), pupilClust);
    printOverlap('Arousal cluster','GLM Pupil', isClustPupil, isGLMpupil);
    printOverlap('Vel cluster A','GLM Vel(+)', isClustVelP, isGLMvelP);
    printOverlap('Vel cluster B','GLM Vel(-)', isClustVelM, isGLMvelM);
end

% ---- map GLM winning feature to group index isn't needed; we use sets ----

function printOverlap(nameA, nameB, A, B)
    both=sum(A&B); onlyA=sum(A&~B); onlyB=sum(~A&B);
    j = both / max(sum(A|B),1);
    fprintf('%s vs %s:  both=%d, %s-only=%d, %s-only=%d, Jaccard=%.2f\n', ...
        nameA, nameB, both, nameA, onlyA, nameB, onlyB, j);
end

% Two-circle Venn with counts. A,B are logical vectors over the same cells.
function drawVenn(A, B, labA, labB, col)
    nBoth=sum(A&B); nA=sum(A&~B); nB=sum(~A&B);
    axis equal off; hold on;
    r=1; d=1.1;                       % circle radius and center separation
    th=linspace(0,2*pi,200);
    cAx=-d/2; cBx=d/2;
    fill(cAx+r*cos(th), r*sin(th), col, 'FaceAlpha',0.35,'EdgeColor',col,'LineWidth',1.5);
    fill(cBx+r*cos(th), r*sin(th), [0.3 0.3 0.3],'FaceAlpha',0.30,'EdgeColor',[0.3 0.3 0.3],'LineWidth',1.5);
    text(cAx-0.4, 0, num2str(nA),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
    text(0, 0, num2str(nBoth),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
    text(cBx+0.4, 0, num2str(nB),'HorizontalAlignment','center','FontSize',12,'FontWeight','bold');
    text(cAx, -1.25, labA,'HorizontalAlignment','center','FontSize',8,'Color',col);
    text(cBx, 1.25, labB,'HorizontalAlignment','center','FontSize',8,'Color',[0.3 0.3 0.3]);
    xlim([-1.8 1.8]); ylim([-1.7 1.7]);
end

function v = nan2neg(x), if isnan(x), v=-Inf; else, v=x; end, end

function s = signLocal(x, y)
    n=min(length(x),length(y)); r=corr(x(1:n),y(1:n),'rows','complete');
    if isnan(r), s=1; else, s=sign(r); if s==0, s=1; end, end
end
