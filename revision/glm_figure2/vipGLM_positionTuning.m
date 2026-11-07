function vipGLM_positionTuning()
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_POSITIONTUNING  Stacked position-tuning figure, grouped by GLM winning
% feature, built from a SAVED cohort cellTable + the session files. Does NOT
% re-run the GLM.
%
% INPUT  STEP 1: COHORT RESULTS file (ONE .mat); STEP 2: SESSION DATA files (one or more).
% OUTPUT a stacked position-tuning figure.
%
%   For each winning-feature group (Vel+, Vel-, Accel, Licking, Reward,
%   Position, Pupil), and each context (Familiar | Novel), shows:
%       top:    mean +/- SEM z-scored dF/F vs track position (pooled cells)
%       bottom: the group's DEFINING feature, mean +/- SEM across trials vs
%               track position, with reward-zone lines (green dashed).
%
%   Inputs (prompted):
%     1) saved cohort results .mat (res.cellTable)
%     2) the session .mat files (multi-select) for dF/F + behaviour
%
%   Cells grouped by their FAMILIAR winning feature (matching the cohort).

    nBins = 40;
    uiwait(msgbox({'STEP 1 of 2:','','Select the SAVED COHORT RESULTS file', ...
        '(cohortGLM_*.mat with res.cellTable).'},'Step 1','help','modal'));
        uiwait(msgbox({'Select the COHORT RESULTS file.','','ONE file: the cohortGLM_*.mat saved by vipGLM_cohort (has res.cellTable).','','NOT a session/m8xxx file.'},'Cohort results (ONE file)','help'));
    [rf,rp]=uigetfile('*.mat','COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes);
    if isequal(rf,0), return; end
    R=load(fullfile(rp,rf));
    if ~isfield(R,'res')||~isfield(R.res,'cellTable'), error('No res.cellTable.'); end
    T=R.res.cellTable; T=T(T.context=="Familiar" & T.ok,:);
    hasVelSign = any(strcmpi(T.Properties.VariableNames,'velSign'));

    uiwait(msgbox({'STEP 2 of 2:','','Select the SESSION .mat files (multi-select).', ...
        'These provide dF/F + behaviour for the position tuning.'},'Step 2','help','modal'));
        uiwait(msgbox({'Select the SESSION DATA files.','','ONE OR MORE m8xxx session .mat files (hold Cmd/Shift for several).'},'Session data (one or more)','help'));
    [sf,sp]=uigetfile('*.mat','SESSION DATA files (ONE OR MORE m8xxx .mat; multi-select)',dfSess,'MultiSelect','on');
    if isequal(sf,0), return; end
    if ischar(sf), sf={sf}; end

    nG=7; catNames={'Vel (+)','Vel (-)','Accel','Licking','Reward','Position','Pupil'};
    % pooled per-cell position-tuning of dF/F (one row per cell), per group/context
    tunFam=cell(1,nG); tunNov=cell(1,nG);
    % per-session feature curves (one row per session), per group/context
    featFam=cell(1,nG); featNov=cell(1,nG);
    rewFamAll={}; rewNovAll={};
    ctr=((1:nBins)-0.5)/nBins;

    for si=1:numel(sf)
        fn=sf{si}; S=load(fullfile(sp,fn));
        if ~isfield(S,'sData'), warning('skip %s',fn); continue; end
        sData=S.sData;
        if ~(isfield(sData.imdata,'roiSignals')&&numel(sData.imdata.roiSignals)>=2)
            warning('skip %s (no roiSignals(2))',fn); continue;
        end
        dff=sData.imdata.roiSignals(2).dff; nCells=size(dff,1); nF=size(dff,2);
        ti=double(sData.behavior.signals.trialIndexesDs); ti=ti(:); N=min(nF,numel(ti)); ti=ti(1:N);
        pos=double(sData.behavior.signals.corridorPositionDs(1:N)); pos=pos(:);
        vel=double(sData.behavior.signals.velocityDs(1:N)); vel=vel(:);
        acc=[0;diff(vel)];
        lick=colvec(getSig(sData,'lickIfreqDs'),N);
        rew =rewardSignal(sData,N);
        pup =pupilTrace(sData,N);
        ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2);
        ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
        vtr=ti>0&ti<=length(ctx); cctx=nan(N,1); cctx(vtr)=ctx(ti(vtr)); cctx=fillmissing(cctx,'nearest');
        fam_id=ctx(1); uctx=unique(ctx); nov_id=uctx(uctx~=fam_id);
        famMask=(cctx==fam_id); novMask=ismember(cctx,nov_id);

        pmin=min(pos); pmax=max(pos); if pmax<=pmin, continue; end
        posN=(pos-pmin)/(pmax-pmin); edges=linspace(0,1,nBins+1);
        binIdx=discretize(posN,edges); binIdx(isnan(binIdx))=1;

        featSig={vel,vel,acc,lick,rew,pos,pup};

        % feature curves per group (same defining signal regardless of cells present)
        for g=1:nG
            featFam{g}=[featFam{g}; featBandMean(featSig{g},binIdx,ti,famMask,nBins)]; %#ok<AGROW>
            featNov{g}=[featNov{g}; featBandMean(featSig{g},binIdx,ti,novMask,nBins)]; %#ok<AGROW>
        end
        rewFamAll{end+1}=rewardPositions(rew,posN,famMask); %#ok<AGROW>
        rewNovAll{end+1}=rewardPositions(rew,posN,novMask); %#ok<AGROW>

        % per-cell dF/F position tuning, grouped by familiar winning feature
        sessRows=T(strcmpi(string(T.session),string(fn)),:);
        for r=1:height(sessRows)
            cid=sessRows.cellID(r); if cid<1||cid>nCells, continue; end
            y=zscoreSafe(double(dff(cid,1:N))');
            g=labelToGroup(sessRows(r,:), hasVelSign, vel, y, famMask);
            if g<1, continue; end
            tunFam{g}=[tunFam{g}; binMean(y,binIdx,famMask,nBins)]; %#ok<AGROW>
            tunNov{g}=[tunNov{g}; binMean(y,binIdx,novMask,nBins)]; %#ok<AGROW>
        end
    end

    rewFam=meanRewLoc(rewFamAll); rewNov=meanRewLoc(rewNovAll);
    c_fam=[0.2 0.5 0.8]; c_nov=[0.8 0.3 0.3]; c_feat=[0.35 0.35 0.35];
    figure('Color','w','Position',[60 40 820 1250]);
    sgtitle('Cohort position tuning (pooled cells) + mean feature; green = reward zones', ...
        'FontWeight','bold','FontSize',12);
    nRow=2*nG;
    for gi=1:nG
        rTop=(gi-1)*4;
        subplot(nRow,2,rTop+1); hold on; drawBand(ctr,tunFam{gi},c_fam);
        ylabel({catNames{gi},'\DeltaF/F (z)'},'FontWeight','bold','FontSize',8);
        if gi==1, title('Familiar'); end
        text(0.02,0.85,sprintf('n=%d',size(tunFam{gi},1)),'Units','normalized','FontSize',8);
        set(gca,'TickDir','out','XTickLabel',[]); box off; xlim([0 1]);
        subplot(nRow,2,rTop+2); hold on; drawBand(ctr,tunNov{gi},c_nov);
        if gi==1, title('Novel'); end
        text(0.02,0.85,sprintf('n=%d',size(tunNov{gi},1)),'Units','normalized','FontSize',8);
        set(gca,'TickDir','out','XTickLabel',[]); box off; xlim([0 1]);
        subplot(nRow,2,rTop+3); hold on; drawBand(ctr,featFam{gi},c_feat); addRewardLines(rewFam);
        ylabel('feature','FontSize',8); if gi==nG, xlabel('Track position (norm.)'); end
        set(gca,'TickDir','out'); box off; xlim([0 1]);
        subplot(nRow,2,rTop+4); hold on; drawBand(ctr,featNov{gi},c_feat); addRewardLines(rewNov);
        if gi==nG, xlabel('Track position (norm.)'); end
        set(gca,'TickDir','out'); box off; xlim([0 1]);
    end
end

function g = labelToGroup(row, hasVelSign, vel, y, famMask)
    g=0; pf=string(row.primaryFeature);
    switch pf
        case "Velocity"
            vs=NaN;
            if hasVelSign && ~isnan(row.velSign), vs=sign(row.velSign); end
            if isnan(vs), vs=signLocal(vel(famMask), y(famMask)); end
            g=(vs>=0)*1 + (vs<0)*2;
        case "Acceleration", g=3;
        case "Licking", g=4;
        case "Reward", g=5;
        case "Position", g=6;
        case "Pupil", g=7;
    end
end

function m = binMean(y, binIdx, mask, nBins)
    m=nan(1,nBins);
    for b=1:nBins, sel=mask&(binIdx==b); if any(sel), m(b)=mean(y(sel),'omitnan'); end, end
end

function m = featBandMean(sig, binIdx, trialIdx, mask, nBins)
    m=nan(1,nBins); if isempty(sig), return; end
    trials=unique(trialIdx(mask & trialIdx>0)); if isempty(trials), return; end
    pt=nan(numel(trials),nBins);
    for tIx=1:numel(trials)
        tsel=mask&(trialIdx==trials(tIx));
        for b=1:nBins, s=tsel&(binIdx==b); if any(s), pt(tIx,b)=mean(sig(s),'omitnan'); end, end
    end
    m=mean(pt,1,'omitnan');
end

function drawBand(x, M, col)
    if isempty(M)||all(isnan(M(:))), text(0.4,0.5,'no cells','Units','normalized','Color',[.5 .5 .5]); return; end
    mu=mean(M,1,'omitnan'); n=sum(~isnan(M),1); sem=std(M,0,1,'omitnan')./sqrt(max(n,1));
    good=~isnan(mu); xx=x(good); mu=mu(good); sem=sem(good);
    fill([xx fliplr(xx)],[mu+sem fliplr(mu-sem)],col,'FaceAlpha',0.25,'EdgeColor','none');
    plot(xx,mu,'Color',col,'LineWidth',1.5);
end

function addRewardLines(locs)
    if isempty(locs), return; end
    yl=ylim; for k=1:numel(locs), plot([locs(k) locs(k)],yl,'--','Color',[0.1 0.6 0.2],'LineWidth',1); end
    ylim(yl);
end

function locs = rewardPositions(rew, posN, mask)
    locs=[]; if isempty(rew), return; end
    ev=mask&(rew(:)>0); if ~any(ev), return; end
    p=sort(posN(ev)); locs=p(1); last=p(1);
    for i=2:numel(p), if p(i)-last>0.03, locs(end+1)=p(i); end, last=p(i); end %#ok<AGROW>
end

function loc = meanRewLoc(rewCells)
    allp=[]; for i=1:numel(rewCells), allp=[allp, rewCells{i}(:)']; end %#ok<AGROW>
    if isempty(allp), loc=[]; return; end
    p=sort(allp); loc=p(1); last=p(1);
    for i=2:numel(p), if p(i)-last>0.05, loc(end+1)=p(i); end, last=p(i); end %#ok<AGROW>
end

function rew = rewardSignal(sData, N)
    rew=zeros(N,1); vt=getSig(sData,'rewardDs');
    if ~isempty(vt)&&any(double(vt)>0), vt=double(vt); rew(1:min(N,numel(vt)))=vt(1:min(N,numel(vt)))>0; return; end
    if isfield(sData,'daqdata')&&isfield(sData.daqdata,'frameIndex')&&isfield(sData.daqdata,'waterValve')
        vr=double(sData.daqdata.waterValve>0.5); ft=find(diff(sData.daqdata.frameIndex)==1)+1; lim=min(N,length(ft)-1);
        for i=1:lim, if any(vr(ft(i):ft(i+1)-1)), rew(i)=1; end, end
    end
end

function v = colvec(x, N)
    if isempty(x), v=zeros(N,1); return; end
    x=double(x(:)); if numel(x)>=N, v=x(1:N); else, v=[x;zeros(N-numel(x),1)]; end
end

function s = getSig(sData, name)
    if isfield(sData.behavior.signals,name), s=sData.behavior.signals.(name); else, s=[]; end
end

function p = pupilTrace(sData, N)
    p=nan(N,1);
    try
        if isfield(sData.behavior,'pupildata') && isfield(sData.behavior.pupildata,'diameter')
            pd=double(sData.behavior.pupildata.diameter); pd(isnan(pd))=inf;
            pfi=find(diff(sData.behavior.signals.pupilCamFrameIndexCorr)==1)+1; ifi=find(diff(sData.daqdata.frameIndex)==1);
            pUs=nan(size(sData.behavior.signals.pupilCamFrameIndexCorr)); vl=min(length(pfi),length(pd)); pUs(pfi(1:vl))=pd(1:vl);
            pUs=fillmissing(pUs,'nearest'); pUs(isinf(pUs))=nan; pr=double(pUs(ifi));
            if length(pr)>=N, p=pr(1:N); else, p=[pr(:);nan(N-length(pr),1)]; end
        end
    catch, end
    p=p(:);
end

function z = zscoreSafe(y)
    s=std(y,'omitnan'); if s<1e-9, z=y-mean(y,'omitnan'); else, z=(y-mean(y,'omitnan'))/s; end
end

function s = signLocal(x, y)
    n=min(length(x),length(y)); r=corr(x(1:n),y(1:n),'rows','complete');
    if isnan(r), s=1; else, s=sign(r); if s==0, s=1; end, end
end
