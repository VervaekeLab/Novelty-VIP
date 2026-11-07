function out = vipGLM_clusterValidation(pathName, files)
% VIPGLM_CLUSTERVALIDATION  Quantitative justification (or refutation) of the
% k-means clustering used in the manuscript PCA figure. Computes, on the SAME
% input the clustering uses (z-scored familiar-block response window -> PCA ->
% first 3 PCs), the standard internal-validation diagnostics:
%
%   (1) PCA variance explained (scree)         -- how low-dimensional is it?
%   (2) Elbow: within-cluster SS vs k (1..8)   -- is there a knee at k=3?
%   (3) Silhouette vs k (2..8)                  -- how separated are clusters?
%   (4) Gap statistic vs a uniform null         -- is ANY k>1 justified at all?
%   (5) Cluster-size stability over restarts    -- are the 3 clusters reproducible?
%
% IMPORTANT (read before using in a rebuttal): these are honest diagnostics, not
% advocacy. If silhouette is low (<~0.25) and the gap statistic favours k=1, the
% data are a CONTINUUM and k=3 is NOT defensible -- presenting these plots would
% strengthen the reviewer's objection, not answer it. Only promote this figure to
% the rebuttal if the numbers genuinely support k=3. The function prints a plain-
% language verdict to the console so you can decide.
%
% STANDALONE: needs only the SESSION DATA files (no cohort .mat). All recorded
% cells are clustered.
%
% Call with no args -> interactive (prompts). Headless: pass (path, files).
%
% INPUT  : the SESSION DATA files -- ONE OR MORE sData .mat (multi-select prompt).
% OUTPUT : a 4-panel validation figure + a plain-language verdict in the console.

    rng(1);
    headless = nargin>=2 && ~isempty(files);
    if headless
        sp=pathName; sf=files; if ischar(sf), sf={sf}; end
    else
        % session DATA live in Desktop/VIP-GCaMP-sData
        dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
        dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData');
        if ~exist(dfSess,'dir'), dfSess=userpath; end
        uiwait(msgbox({'Select the SESSION DATA files.','', ...
            'ONE OR MORE m8xxx session .mat files', ...
            '(hold Cmd/Shift to select several at once).'},'Session data (one or more)','help'));
        [sf,sp]=uigetfile('*.mat','SESSION DATA files (ONE OR MORE m8xxx .mat; multi-select)',dfSess,'MultiSelect','on');
        if isequal(sf,0), return; end
        if ischar(sf), sf={sf}; end
    end

    % ---- build the feature matrix the clustering uses (matches pcaOverlay) ----
    % binnedRoisDff is [trials x bins x cells]; reshape to [trials*bins x cells]
    % so each COLUMN is a cell. All recorded cells are included (standalone: no
    % external cell list needed -- the GLM 'ok' filter excluded no cells anyway).
    binnedFiringDff=[];
    for si=1:numel(sf)
        S=load(fullfile(sp,sf{si}));
        if ~isfield(S,'sData'), continue; end
        if ~(isfield(S.sData,'imdata') && isfield(S.sData.imdata,'binnedRoisDff')), continue; end
        B=S.sData.imdata.binnedRoisDff;
        if iscell(B), B=cat(3,B{:}); end
        if ndims(B)~=3, continue; end          % require [trials x bins x cells]
        sz=size(B);
        Bre=reshape(permute(B,[2 1 3]),[sz(1)*sz(2),sz(3)]);   % (trials*bins) x cells
        binnedFiringDff=[binnedFiringDff, Bre]; %#ok<AGROW>
    end
    if isempty(binnedFiringDff), error('No cells found (sessions missing imdata.binnedRoisDff?).'); end

    data=binnedFiringDff(:,:)';
    binIndex=1:size(binnedFiringDff,1); event_time=125*30+1;
    B1_idx=find(binIndex>=30 & binIndex<event_time);
    wd=data(:,B1_idx);
    wz=(wd-mean(wd,2))./std(wd,0,2); wz(isnan(wz))=0;
    [coeff,score,latent]=pca(wz); %#ok<ASGLU>
    nPCs=min(3,size(score,2));
    X=score(:,1:nPCs);
    N=size(X,1);

    % ---- (1) scree / variance explained ----
    ve=100*latent/sum(latent); cve=cumsum(ve);

    % ---- (2) elbow: within-cluster sum of squares ----
    kmax=min(8,N-1); ks=1:kmax; wss=nan(1,kmax);
    for k=ks
        if k==1
            wss(k)=sum(sum((X-mean(X,1)).^2));
        else
            [~,~,sumd]=kmeans(X,k,'Replicates',20,'Display','off','EmptyAction','singleton');
            wss(k)=sum(sumd);
        end
    end

    % ---- (3) silhouette vs k ----
    sil=nan(1,kmax);
    for k=2:kmax
        cl=kmeans(X,k,'Replicates',20,'Display','off','EmptyAction','singleton');
        sil(k)=mean(silhouette(X,cl));
    end

    % ---- (4) gap statistic vs uniform null (the key test: is k>1 justified?) ----
    nRef=50; gap=nan(1,kmax); sk=nan(1,kmax);
    lo=min(X,[],1); hi=max(X,[],1);
    logW=nan(1,kmax); logWref=nan(nRef,kmax);
    for k=1:kmax
        logW(k)=log(localWk(X,k));
        for r=1:nRef
            Xr=lo+(hi-lo).*rand(N,nPCs);
            logWref(r,k)=log(localWk(Xr,k));
        end
        gap(k)=mean(logWref(:,k))-logW(k);
        sk(k)=std(logWref(:,k))*sqrt(1+1/nRef);
    end
    % first k where gap(k) >= gap(k+1)-s(k+1)  (Tibshirani criterion)
    kGap=kmax;
    for k=1:kmax-1
        if gap(k)>=gap(k+1)-sk(k+1), kGap=k; break; end
    end

    % ---- (5) cluster-size stability at k=3 over restarts ----
    nRuns=50; sizes3=zeros(nRuns,3);
    for r=1:nRuns
        cl=kmeans(X,3,'Replicates',1,'Start','plus','Display','off','EmptyAction','singleton');
        s=sort(histcounts(cl,1:4),'descend'); sizes3(r,:)=s;
    end
    sizeMean=mean(sizes3,1); sizeStd=std(sizes3,0,1);

    out=struct('ve',ve,'cve',cve,'ks',ks,'wss',wss,'sil',sil, ...
        'gap',gap,'sk',sk,'kGap',kGap,'sil3',sil(min(3,kmax)), ...
        'sizeMean',sizeMean,'sizeStd',sizeStd,'N',N);

    % ===================== verdict (printed) =====================
    fprintf('\n================ k=3 CLUSTERING VALIDATION (N=%d cells) ================\n',N);
    fprintf('PCA: PC1-3 explain %.1f%% of variance (%.1f / %.1f / %.1f).\n',cve(min(3,end)),ve(1),ve(2),ve(min(3,end)));
    fprintf('Silhouette at k=3: %.3f  ', sil(min(3,kmax)));
    if sil(min(3,kmax))>=0.5,      fprintf('(strong separation)\n');
    elseif sil(min(3,kmax))>=0.25, fprintf('(weak/structure-ish)\n');
    else,                          fprintf('(NO meaningful separation)\n'); end
    [~,kSilBest]=max(sil); fprintf('Best silhouette is at k=%d.\n',kSilBest);
    fprintf('Gap statistic favours k=%d (vs a no-cluster uniform null).\n',kGap);
    fprintf('Cluster sizes at k=3 over %d restarts: %.0f / %.0f / %.0f  (+/- %.0f / %.0f / %.0f)\n',...
        nRuns,sizeMean,sizeStd);
    fprintf('-------------------------------------------------------------------------\n');
    if kGap<=1 || sil(min(3,kmax))<0.25
        fprintf('VERDICT: data look like a CONTINUUM. k=3 is NOT well supported.\n');
        fprintf('         Do NOT present these as a defence of k=3 -- they confirm the\n');
        fprintf('         reviewer. Reframe PCA as a continuum coloured by GLM instead.\n');
    elseif kGap==3 && sil(min(3,kmax))>=0.25
        fprintf('VERDICT: k=3 is genuinely supported (gap picks 3, silhouette ok).\n');
        fprintf('         These plots are a valid rebuttal.\n');
    else
        fprintf('VERDICT: MIXED. Gap picks k=%d, best silhouette k=%d. k=3 is a choice,\n',kGap,kSilBest);
        fprintf('         not a clear optimum -- present honestly with caveats.\n');
    end
    fprintf('=========================================================================\n\n');

    % ===================== figure =====================
    figure('Color','w','Position',[60 60 1100 760],'Name','k=3 clustering validation');
    LW=1.4;
    subplot(2,2,1);
    bar(1:numel(ve),ve,'FaceColor',[0.6 0.6 0.6],'EdgeColor','k'); hold on;
    plot(1:numel(cve),cve,'k-o','LineWidth',LW,'MarkerFaceColor','k');
    xlim([0.3 min(10,numel(ve))+0.7]); xlabel('Principal component'); ylabel('% variance (bar) / cumulative (line)');
    title('(1) PCA scree'); set(gca,'TickDir','out','LineWidth',LW); box off;

    subplot(2,2,2);
    plot(ks,wss,'k-o','LineWidth',LW,'MarkerFaceColor','k'); hold on;
    plot([3 3],ylim,'r--','LineWidth',1);
    xlabel('k (clusters)'); ylabel('Within-cluster SS'); title('(2) Elbow');
    set(gca,'TickDir','out','LineWidth',LW); box off;

    subplot(2,2,3);
    plot(2:kmax,sil(2:kmax),'k-o','LineWidth',LW,'MarkerFaceColor','k'); hold on;
    plot([3 3],ylim,'r--','LineWidth',1); plot(xlim,[0.25 0.25],'-','Color',[.6 .6 .6]);
    xlabel('k (clusters)'); ylabel('Mean silhouette'); title('(3) Silhouette');
    set(gca,'TickDir','out','LineWidth',LW); box off;

    subplot(2,2,4);
    errorbar(ks,gap,sk,'k-o','LineWidth',LW,'MarkerFaceColor','k'); hold on;
    plot([kGap kGap],ylim,'b--','LineWidth',1); plot([3 3],ylim,'r--','LineWidth',1);
    xlabel('k (clusters)'); ylabel('Gap statistic'); title(sprintf('(4) Gap',kGap));
    set(gca,'TickDir','out','LineWidth',LW); box off;

    sgtitle(sprintf('k-means validation  ',...
        cve(min(3,end)),sil(min(3,kmax)),kGap),'FontWeight','bold');
end

function W = localWk(X,k)
    % pooled within-cluster dispersion (sum of within-cluster SS)
    if k==1
        W=sum(sum((X-mean(X,1)).^2)); return;
    end
    [~,~,sumd]=kmeans(X,k,'Replicates',10,'Display','off','EmptyAction','singleton');
    W=sum(sumd);
end
