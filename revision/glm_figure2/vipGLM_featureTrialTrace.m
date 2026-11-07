function vipGLM_featureTrialTrace()
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_FEATURETRIALTRACE  Panel-j-style dF/F vs trial, grouped by GLM feature.
%
% INPUT  STEP 1: COHORT RESULTS file (ONE .mat); STEP 2: SESSION DATA files (one or more).
%
%   Tests whether GLM-identified feature groups (especially Pupil) show the
%   large session-onset / novel-onset transient with slow decay seen in
%   Fig.1 panel j and pupil panel h -- i.e. whether the "pupil/arousal" cells
%   are the onset-transient population.
%
%   Reads:
%     1) a saved cohort results .mat (res.cellTable) for the per-cell GLM
%        winning-feature classification, and
%     2) the original session .mat files for the raw dF/F traces.
%   Cells are matched by session name + cellID. Activity is z-scored per cell
%   (matching the manuscript's z-scored F/F), then averaged across each group's
%   cells and plotted as a dense per-frame trace over the full F->N->F session.
%   Sessions are aligned to novel-block onset and resampled to a common trial
%   axis before averaging (so onset transients line up across sessions).
%
%   Grouping uses each cell's FAMILIAR winning feature (the onset transient
%   lives in the familiar block).

    % --- results file (classification) ---
    uiwait(msgbox({'STEP 1 of 2:', '', ...
        'Select the SAVED COHORT RESULTS file', ...
        '(the cohortGLM_*.mat produced by vipGLM_cohort,', ...
        'which contains res.cellTable with the GLM classification).'}, ...
        'Step 1: results file','help','modal'));
    [rf,rp]=uigetfile('*.mat','STEP 1 of 2 -- select the COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes);
    if isequal(rf,0), return; end
    R=load(fullfile(rp,rf));
    if ~isfield(R,'res')||~isfield(R.res,'cellTable'), error('No res.cellTable in that file.'); end
    T=R.res.cellTable;
    T=T(T.context=="Familiar" & T.ok,:);   % classification from familiar block

    % --- session files (raw dF/F) ---
    uiwait(msgbox({'STEP 2 of 2:', '', ...
        'Now select the SESSION data files', ...
        '(the original sData .mat files that were analyzed -', ...
        'you can select multiple at once with Shift/Cmd-click).', '', ...
        'These provide the raw dF/F traces.'}, ...
        'Step 2: session files','help','modal'));
    [sf,sp]=uigetfile('*.mat','STEP 2 of 2 -- select the SESSION DATA files (ONE OR MORE sData .mat; multi-select)', ...
        dfSess,'MultiSelect','on');
    if isequal(sf,0), return; end
    if ischar(sf), sf={sf}; end

    feats={'Velocity','Acceleration','Licking','Reward','Position','Pupil'};
    nG=7; catNames={'Vel (+)','Vel (-)','Accel','Licking','Reward','Position','Pupil'};
    nTrialGrid=110;        % common trial axis (like panel j range)
    align=cell(1,nG);      % per group: matrix of [resampled-trace] rows (one per cell)
    for g=1:nG, align{g}=[]; end
    novOnsetFrac=[];       % record novel-onset position on the common grid

    for si=1:numel(sf)
        fn=sf{si};
        S=load(fullfile(sp,fn));
        if ~isfield(S,'sData'), warning('skip %s (no sData)',fn); continue; end
        sData=S.sData;
        if ~(isfield(sData.imdata,'roiSignals')&&numel(sData.imdata.roiSignals)>=2)
            warning('skip %s (no roiSignals(2))',fn); continue;
        end
        dff=sData.imdata.roiSignals(2).dff; nCells=size(dff,1); nF=size(dff,2);

        % trial index + context per frame
        ti=double(sData.behavior.signals.trialIndexesDs); ti=ti(:); N=min(nF,numel(ti));
        ti=ti(1:N);
        ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2);
        ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
        vtr=ti>0&ti<=length(ctx); cctx=nan(N,1); cctx(vtr)=ctx(ti(vtr)); cctx=fillmissing(cctx,'nearest');
        fam_id=ctx(1); uctx=unique(ctx); nov_id=uctx(uctx~=fam_id);
        if isempty(nov_id), warning('skip %s (no novel context)',fn); continue; end
        famMask = (cctx==fam_id);   % for velocity-sign computation
        vel = double(sData.behavior.signals.velocityDs(1:N)); vel=vel(:);

        % per-trial mean z-scored dF/F per cell, then we keep trial axis
        trials=unique(ti(ti>0));
        novFirstTrial = min(ti(cctx==nov_id(1) & ti>0));
        % match this session's cells to the classification table
        sessRows = T(strcmpi(string(T.session), string(fn)), :);
        for r=1:height(sessRows)
            cid=sessRows.cellID(r); if cid<1||cid>nCells, continue; end
            g=featToGroup(sessRows(r,:), feats);   % winning feature -> group 1 (Vel), 3..7
            if g<1, continue; end
            y=double(dff(cid,1:N))'; y=zscoreSafe(y);
            % velocity sign-split: group 1 -> 1 (+) or 2 (-), computed in familiar block
            if g==1
                vs = signCorrLocal(vel(famMask), y(famMask));
                if vs<0, g=2; end
            end
            % per-trial mean
            tv=nan(1,numel(trials));
            for k=1:numel(trials)
                sel=ti==trials(k); if any(sel), tv(k)=mean(y(sel),'omitnan'); end
            end
            % build a trial axis centered so novel onset is at a fixed fraction
            tn = double(trials(:)') ;                  % trial numbers present
            % resample onto common grid 1..nTrialGrid with novel onset at ~0.27
            grid = resampleToGrid(tn, tv, novFirstTrial, nTrialGrid);
            align{g}=[align{g}; grid];
        end
        novOnsetFrac(end+1)=0.27; %#ok<AGROW>
    end

    % --- plot ---
    colors=[0.85 0.33 0.10; 0.93 0.69 0.13; 0.00 0.45 0.74; ...
            0.49 0.18 0.56; 0.47 0.67 0.19; 0.30 0.75 0.93; 0.64 0.08 0.18];
    novStart=round(0.27*nTrialGrid); novEnd=round(0.82*nTrialGrid);
    figure('Color','w','Position',[80 150 1150 520]); hold on;

    % shading for novel block
    h=[]; labels={}; allv=[];
    for g=1:nG
        M=align{g}; if isempty(M), continue; end
        mu=mean(M,1,'omitnan'); sem=std(M,0,1,'omitnan')./sqrt(max(sum(~isnan(M),1),1));
        x=1:numel(mu); good=~isnan(mu);
        fill([x(good) fliplr(x(good))],[mu(good)+sem(good) fliplr(mu(good)-sem(good))], ...
            colors(g,:),'FaceAlpha',0.15,'EdgeColor','none');
        hh=plot(x(good),mu(good),'Color',colors(g,:),'LineWidth',1.5);
        h(end+1)=hh; labels{end+1}=sprintf('%s (n=%d)',catNames{g},size(M,1)); %#ok<AGROW>
        allv=[allv mu(good)]; %#ok<AGROW>
    end
    if ~isempty(allv)
        ylo=min(allv)*1.1; yhi=max(allv)*1.1;
        pa=patch([novStart novEnd novEnd novStart],[ylo ylo yhi yhi],[0.9 0.9 0.9],'EdgeColor','none');
        uistack(pa,'bottom'); ylim([ylo yhi]);
    end
    xlabel('Trial number (aligned to novel onset)'); ylabel('Z-scored \DeltaF/F (mean \pm SEM)');
    title('VIP dF/F vs trial, grouped by GLM winning feature (familiar-block classification)','FontWeight','bold');
    if ~isempty(h), legend(h,labels,'Location','northeast','FontSize',8); end
    box off; set(gca,'TickDir','out'); xlim([1 nTrialGrid]);
end

% winning feature -> group index (velocity sign-split into 1/2 using stored dR2 sign? 
% the table doesn't store sign, so velocity maps to group 1 (+) by default; the
% sign split is not needed for the transient question).
function g = featToGroup(row, feats)
    g=0; pf=string(row.primaryFeature);
    switch pf
        case "Velocity", g=1;
        case "Acceleration", g=3;
        case "Licking", g=4;
        case "Reward", g=5;
        case "Position", g=6;
        case "Pupil", g=7;
        otherwise, g=0;
    end
end

function grid = resampleToGrid(tn, tv, novFirstTrial, nGrid)
    % Map trial numbers so that novFirstTrial lands at 0.27*nGrid; linear before
    % and after. Then interpolate tv onto 1..nGrid.
    grid=nan(1,nGrid);
    good=~isnan(tv); if nnz(good)<3, return; end
    tn=tn(good); tv=tv(good);
    % normalized trial coordinate: 0 at first trial, with novel onset at 0.27
    t0=min(tn); tNov=novFirstTrial; tEnd=max(tn);
    if tNov<=t0 || tEnd<=tNov, % degenerate; just stretch linearly
        u=(tn-t0)/max(tEnd-t0,1);
    else
        u=zeros(size(tn));
        pre=tn<tNov; u(pre)=0.27*(tn(pre)-t0)/(tNov-t0);
        post=~pre;  u(post)=0.27+0.73*(tn(post)-tNov)/(tEnd-tNov);
    end
    xq=linspace(0,1,nGrid);
    grid=interp1(u,tv,xq,'linear',NaN);
end

function z = zscoreSafe(y)
    s=std(y,'omitnan'); if s<1e-9, z=y-mean(y,'omitnan'); else, z=(y-mean(y,'omitnan'))/s; end
end

function s = signCorrLocal(x, y)
    n=min(length(x),length(y));
    r=corr(x(1:n),y(1:n),'rows','complete');
    if isnan(r), s=1; else, s=sign(r); if s==0, s=1; end, end
end
