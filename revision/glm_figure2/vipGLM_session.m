function vipGLM_session()
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_SESSION  Population GLM summary for ONE VIP session (10-30 cells).
%
% AIM: a four-panel summary (significant modulations, primary feature,
% multiplexing distribution, effect-size dR2), split Familiar vs Novel with
% velocity/pupil sign-splitting, for a single session. Useful for per-session QC.
% Uses the same engine (vipGLM_fitCell) as the cohort, so results are identical.
%
% INPUT:  ONE SESSION DATA file (sData .mat), via prompt.
% OUTPUT: one summary figure. Requires vipGLM_fitCell.m + vipGLM_defaultFolder.m.

    userHome = getenv('HOME');
    folderPath = fullfile(userHome,'Desktop','VIP-GCaMP-sData');
    if ~exist(folderPath,'dir'), folderPath = fullfile(userHome,'Desktop'); end
        uiwait(msgbox({'Select ONE SESSION DATA file.','','A single m8xxx session .mat file.'},'Session data (one file)','help'));
    [fileName,pathName] = uigetfile('*.mat','ONE SESSION DATA file (m8xxx .mat) to summarize',dfSess);
    if isequal(fileName,0), return; end

    fprintf('Loading %s...\n', fileName);
    S = load(fullfile(pathName,fileName));
    sData = S.sData;
    if ~(isfield(sData.imdata,'roiSignals') && numel(sData.imdata.roiSignals)>=2)
        error('Could not find VIP neural data in sData.imdata.roiSignals(2).dff');
    end
    dff_all = sData.imdata.roiSignals(2).dff;

    uiFig = uifigure('Name',['Batch GLM: ' fileName],'Position',[100 100 380 430]);
    g = uigridlayout(uiFig,[10 2]); g.RowHeight=repmat({'fit'},1,10); g.ColumnWidth={150,'1x'};
    uilabel(g,'Text','Ridge \lambda:');         lamF  = uieditfield(g,'numeric','Value',0.5);
    uilabel(g,'Text','Lambda mode:');            lamDD = uidropdown(g,'Items',{'Fixed lambda','Nested CV lambda'},'Value','Fixed lambda');
    uilabel(g,'Text','Significance:');           sigDD = uidropdown(g,'Items',{'Permutation (rigorous)','t-test (fast)'},'Value','Permutation (rigorous)');
    uilabel(g,'Text','# permutations:');         permF = uieditfield(g,'numeric','Value',200,'Limits',[20 2000],'RoundFractionalValues','on');
    uilabel(g,'Text','Lag min (s):');            lagMinF = uieditfield(g,'numeric','Value',-3,'Limits',[-10 10]);
    uilabel(g,'Text','Lag max (s):');            lagMaxF = uieditfield(g,'numeric','Value',3,'Limits',[-10 10]);
    uilabel(g,'Text','# lag points:');           nLagF   = uieditfield(g,'numeric','Value',5,'Limits',[1 21],'RoundFractionalValues','on');
    uilabel(g,'Text','Use parfor:');             parDD = uidropdown(g,'Items',{'Yes','No'},'Value','Yes');
    uilabel(g,'Text','Save results .mat:');      saveDD= uidropdown(g,'Items',{'Yes','No'},'Value','No');
    runB = uibutton(g,'Text','Run Population GLM','ButtonPushedFcn', @(s,e) runBatchGLM( ...
        sData, dff_all, lamF.Value, round(permF.Value), strcmp(parDD.Value,'Yes'), strcmp(saveDD.Value,'Yes'), ...
        ternary(contains(lamDD.Value,'Nested'),'nested','fixed'), ternary(contains(sigDD.Value,'t-test'),'ttest','perm'), ...
        lagMinF.Value, lagMaxF.Value, round(nLagF.Value), uiFig, fileName));
    runB.Layout.Row=10; runB.Layout.Column=[1 2];
end

function v = ternary(cond,a,b), if cond, v=a; else, v=b; end, end

% =====================================================================
function runBatchGLM(sData, dff_all, lambdaVal, nPerms, usePar, doSave, lambdaMode, sigTest, lagMin, lagMax, nLags, uiFig, fileName)
    nCells = size(dff_all,1);
    d = uiprogressdlg(uiFig,'Title','Population GLM','Message','Fitting cells...','Value',0.02);
    cleanup = onCleanup(@() closeIfValid(d)); %#ok<NASGU>

    % velocity & pupil traces (only needed to determine +/- sign of modulation)
    [velTrace, pupilTrace, ~, N, famMask, novMask] = quickBehavior(sData, size(dff_all,2));

    % ---- fit every cell in BOTH contexts via the hardened engine ----
    famRes = cell(nCells,1); novRes = cell(nCells,1);
    seeds = 1:nCells;
    runPar = usePar && ~isempty(ver('parallel'));

    if runPar
        dq = parallel.pool.DataQueue; progressTick(d,nCells,true);
        afterEach(dq,@(~) progressTick(d,nCells,false));
        d.Message = sprintf('Fitting %d cells x 2 contexts on pool...', nCells);
        parfor c=1:nCells
            rng(seeds(c));
            famRes{c} = vipGLM_fitCell(sData, dff_all, c, 'Familiar Context', lambdaVal, nPerms, lambdaMode, sigTest, lagMin,lagMax,nLags);
            novRes{c} = vipGLM_fitCell(sData, dff_all, c, 'Novel Context',    lambdaVal, nPerms, lambdaMode, sigTest, lagMin,lagMax,nLags);
            send(dq,c);
        end
    else
        for c=1:nCells
            rng(seeds(c));
            famRes{c} = vipGLM_fitCell(sData, dff_all, c, 'Familiar Context', lambdaVal, nPerms, lambdaMode, sigTest, lagMin,lagMax,nLags);
            novRes{c} = vipGLM_fitCell(sData, dff_all, c, 'Novel Context',    lambdaVal, nPerms, lambdaMode, sigTest, lagMin,lagMax,nLags);
            d.Value = 0.02 + 0.93*(c/nCells); d.Message = sprintf('Fitting cell %d/%d...', c, nCells);
        end
    end

    % ---- aggregate into the 8 sign-split categories ----
    n_cats = 7;
    cat_names = {'Vel (+)','Vel (-)','Accel','Licking','Reward','Position','Pupil'};
    sig_counts_fam=zeros(1,n_cats); prim_counts_fam=zeros(1,n_cats); dR2_vals_fam=cell(1,n_cats);
    sig_counts_nov=zeros(1,n_cats); prim_counts_nov=zeros(1,n_cats); dR2_vals_nov=cell(1,n_cats);
    num_feats_fam=zeros(nCells,1); num_feats_nov=zeros(nCells,1);
    primCat_fam=zeros(nCells,1); primCat_nov=zeros(nCells,1);   % winning category 1-7 (0=none)

    for c=1:nCells
        % sign of velocity / pupil modulation, computed WITHIN each context
        ydff = double(dff_all(c,1:N))';
        velSign_f   = signCorr(velTrace(famMask),   ydff(famMask));
        velSign_n   = signCorr(velTrace(novMask),   ydff(novMask));
        pupilSign_f = signCorr(pupilTrace(famMask), ydff(famMask));
        pupilSign_n = signCorr(pupilTrace(novMask), ydff(novMask));

        of = famRes{c};
        if ~isempty(of) && of.ok
            [sig_counts_fam,prim_counts_fam,dR2_vals_fam,nf] = ...
                tally(of, velSign_f, pupilSign_f, of.hasPupil, sig_counts_fam, prim_counts_fam, dR2_vals_fam);
            num_feats_fam(c)=nf;
            primCat_fam(c)=primaryCat(of, velSign_f);
        end
        on = novRes{c};
        if ~isempty(on) && on.ok
            [sig_counts_nov,prim_counts_nov,dR2_vals_nov,nf] = ...
                tally(on, velSign_n, pupilSign_n, on.hasPupil, sig_counts_nov, prim_counts_nov, dR2_vals_nov);
            num_feats_nov(c)=nf;
            primCat_nov(c)=primaryCat(on, velSign_n);
        end
    end

    % ---- summary stats ----
    mean_dR2_fam=zeros(1,n_cats); err_dR2_fam=zeros(1,n_cats);
    mean_dR2_nov=zeros(1,n_cats); err_dR2_nov=zeros(1,n_cats);
    for i=1:n_cats
        if ~isempty(dR2_vals_fam{i})
            mean_dR2_fam(i)=mean(dR2_vals_fam{i},'omitnan'); err_dR2_fam(i)=std(dR2_vals_fam{i},'omitnan')/sqrt(numel(dR2_vals_fam{i}));
        end
        if ~isempty(dR2_vals_nov{i})
            mean_dR2_nov(i)=mean(dR2_vals_nov{i},'omitnan'); err_dR2_nov(i)=std(dR2_vals_nov{i},'omitnan')/sqrt(numel(dR2_vals_nov{i}));
        end
    end
    dist_fam = histcounts(num_feats_fam, -0.5:6.5);
    dist_nov = histcounts(num_feats_nov, -0.5:6.5);

    % ---- PLOTS (kept from the original) ----
    figure('Name',['Population GLM: ' fileName],'Position',[100 100 1400 800],'Color','w');
    if strcmpi(lambdaMode,'nested'), lamStr='\lambda=CV-tuned'; else, lamStr=sprintf('\\lambda=%.3g',lambdaVal); end
    sgtitle(sprintf('Population Encoding (%d VIP cells) | %s | %d perms\n%s', ...
        nCells, lamStr, nPerms, strrep(fileName,'_','\_')),'FontWeight','bold','FontSize',14);
    c_fam=[0.2 0.5 0.8]; c_nov=[0.8 0.3 0.3];

    subplot(2,2,1); hold on;
    b1=bar(1:n_cats,[sig_counts_fam' sig_counts_nov'],'grouped','EdgeColor','k');
    b1(1).FaceColor=c_fam; b1(2).FaceColor=c_nov;
    set(gca,'XTick',1:n_cats,'XTickLabel',cat_names,'XTickLabelRotation',30,'TickDir','out');
    ylabel('Number of VIP cells'); title('1. Significant modulations (p<0.05, \DeltaR^2>0.01)'); box off;
    legend({'Familiar','Novel'},'Location','best');

    subplot(2,2,2); hold on;
    b2=bar(1:n_cats,[prim_counts_fam' prim_counts_nov'],'grouped','EdgeColor','k');
    b2(1).FaceColor=c_fam; b2(2).FaceColor=c_nov;
    set(gca,'XTick',1:n_cats,'XTickLabel',cat_names,'XTickLabelRotation',30,'TickDir','out');
    ylabel('Number of VIP cells'); title('2. Primary (winning) feature per cell'); box off;

    subplot(2,2,3); hold on;
    b3=bar(0:6,[dist_fam' dist_nov'],'grouped','EdgeColor','k');
    b3(1).FaceColor=c_fam; b3(2).FaceColor=c_nov;
    xlabel('Number of significant features per cell'); ylabel('Number of VIP cells');
    title('3. Multiplexing distribution (0 = unresponsive or shared-only)'); set(gca,'TickDir','out'); box off;

    subplot(2,2,4); hold on;
    b4=bar(1:n_cats,[mean_dR2_fam' mean_dR2_nov'],'grouped','EdgeColor','k');
    b4(1).FaceColor=c_fam; b4(2).FaceColor=c_nov;
    for i=1:2
        if i==1, ed=err_dR2_fam; else, ed=err_dR2_nov; end
        errorbar(b4(i).XEndPoints,b4(i).YData,ed,'k','LineStyle','none','LineWidth',1);
    end
    set(gca,'XTick',1:n_cats,'XTickLabel',cat_names,'XTickLabelRotation',30,'TickDir','out');
    ylabel('Mean unique \DeltaR^2 (\pm SEM)'); title('4. Effect size (significant cells only)'); box off;

    % ---- grouped position-tuning figure: cells grouped by winning feature ----
    plotGroupedPositionTuning(dff_all, sData, famMask, novMask, N, primCat_fam, primCat_nov, cat_names, fileName);

    if doSave
        outFile=sprintf('popGLM_%s_%s.mat', erase(fileName,'.mat'), datestr(now,'yyyymmdd_HHMMSS'));
        res=struct('famRes',{famRes},'novRes',{novRes},'cat_names',{cat_names}, ...
            'sig_fam',sig_counts_fam,'sig_nov',sig_counts_nov,'prim_fam',prim_counts_fam,'prim_nov',prim_counts_nov, ...
            'dist_fam',dist_fam,'dist_nov',dist_nov,'meta',struct('lambda',lambdaVal,'lambdaMode',lambdaMode,'nPerms',nPerms,'lagMin',lagMin,'lagMax',lagMax,'nLags',nLags,'file',fileName)); %#ok<NASGU>
        save(outFile,'res'); fprintf('Saved %s\n',outFile);
    end
    close(d);
end

% =====================================================================
% Map a vipGLM_fitCell result (6 features) into the 8 sign-split categories.
% feature order from vipGLM_fitCell: 1 Vel, 2 Accel, 3 Lick, 4 Reward, 5 Position, 6 Pupil
function [s_counts, p_counts, dr2_vals, nFeat] = tally(o, velSign, pupilSign, hasPupil, s_counts, p_counts, dr2_vals)
    sig = o.sig; dR2 = o.dR2_unique;
    nFeat = sum(sig(~isnan(dR2)));   % count significant features (pupil NaN if absent)

    % primary = largest unique dR2 among significant features
    primIdx = 0; valid = find(sig & ~isnan(dR2));
    if ~isempty(valid), [~,m]=max(dR2(valid)); primIdx = valid(m); end

    % Velocity (feature 1) -> category 1 (+) or 2 (-)
    if sig(1)
        idx = (velSign>=0)*1 + (velSign<0)*2;
        s_counts(idx)=s_counts(idx)+1; dr2_vals{idx}=[dr2_vals{idx}, dR2(1)];
        if primIdx==1, p_counts(idx)=p_counts(idx)+1; end
    end
    % Accel(2)->3, Lick(3)->4, Reward(4)->5, Position(5)->6
    map = [3 4 5 6];
    for f=2:5
        if sig(f)
            idx=map(f-1);
            s_counts(idx)=s_counts(idx)+1; dr2_vals{idx}=[dr2_vals{idx}, dR2(f)];
            if primIdx==f, p_counts(idx)=p_counts(idx)+1; end
        end
    end
    % Pupil (feature 6) -> category 7 (single, no sign split)
    if hasPupil && sig(6) && ~isnan(dR2(6))
        idx = 7;
        s_counts(idx)=s_counts(idx)+1; dr2_vals{idx}=[dr2_vals{idx}, dR2(6)];
        if primIdx==6, p_counts(idx)=p_counts(idx)+1; end
    end
end


% Winning (primary) category 1-7 for a cell, matching the tally scheme:
% 1 Vel(+) 2 Vel(-) 3 Accel 4 Licking 5 Reward 6 Position 7 Pupil. 0 = none.
function cat = primaryCat(o, velSign)
    cat = 0;
    sig=o.sig; dR2=o.dR2_unique;
    valid=find(sig & ~isnan(dR2));
    if isempty(valid), return; end
    [~,m]=max(dR2(valid)); pf=valid(m);
    switch pf
        case 1, cat=(velSign>=0)*1 + (velSign<0)*2;   % velocity sign-split
        case 2, cat=3;   % accel
        case 3, cat=4;   % licking
        case 4, cat=5;   % reward
        case 5, cat=6;   % position
        case 6, cat=7;   % pupil
    end
end

% =====================================================================
% Group cells by winning feature and plot mean +/- SEM position-tuning
% (dF/F averaged across trials vs track position), Familiar | Novel columns.
function plotGroupedPositionTuning(dff_all, sData, famMask, novMask, N, primCat_fam, primCat_nov, cat_names, fileName)
    nBins = 40;
    pos = double(sData.behavior.signals.corridorPositionDs(1:N)); pos = pos(:);
    pmin=min(pos); pmax=max(pos); if pmax<=pmin, return; end
    posN=(pos-pmin)/(pmax-pmin); edges=linspace(0,1,nBins+1); ctr=(edges(1:end-1)+edges(2:end))/2;
    binIdx=discretize(posN,edges); binIdx(isnan(binIdx))=1;
    trialIdx = double(sData.behavior.signals.trialIndexesDs(1:N)); trialIdx=trialIdx(:);

    nCells=size(dff_all,1);
    % per-cell position-binned mean dF/F within each context
    tun_fam=nan(nCells,nBins); tun_nov=nan(nCells,nBins);
    for c=1:nCells
        y=double(dff_all(c,1:N))'; y=zscore_safe(y);
        tun_fam(c,:)=binMean(y, binIdx, famMask, nBins);
        tun_nov(c,:)=binMean(y, binIdx, novMask, nBins);
    end

    % behavioral feature traces (raw signals, per category 1-7)
    vel = colvec(sData.behavior.signals.velocityDs, N);
    acc = [0; diff(vel)];
    lick = colvec(getSig(sData,'lickIfreqDs'), N);
    rew  = rewardSignal(sData, N);
    pup  = pupilTraceLocal(sData, N);
    % map each category to its defining behavioral signal (position group -> [] = occupancy)
    featSig = {vel, vel, acc, lick, rew, pos, pup};   % cats 1..7 (vel+/-, accel, lick, reward, position, pupil)

    % reward locations: normalized track positions where reward was delivered
    rewLocFam = rewardPositions(rew, posN, famMask);
    rewLocNov = rewardPositions(rew, posN, novMask);

    groups=1:7; nG=numel(groups);
    figure('Name',['Position tuning by winning feature: ' fileName], ...
        'Color','w','Position',[80 60 760 1200]);
    sgtitle(sprintf('Position tuning (top) + defining feature (bottom), grouped by winning feature\n%s', ...
        strrep(fileName,'_','\_')),'FontWeight','bold','FontSize',12);
    c_fam=[0.2 0.5 0.8]; c_nov=[0.8 0.3 0.3]; c_feat=[0.35 0.35 0.35];
    nRow = 2*nG;   % each group: dF/F row then feature row

    for gi=1:nG
        g=groups(gi);
        rTop=(gi-1)*2*2;   % subplot index base for this group (2 cols x 2 subrows)
        % --- dF/F top, Familiar (col 1) ---
        subplot(nRow,2,rTop+1); hold on;
        cellsF=find(primCat_fam==g);
        drawBand(ctr, tun_fam(cellsF,:), c_fam);
        ylabel({cat_names{g},'\DeltaF/F (z)'},'FontWeight','bold','FontSize',8);
        if gi==1, title('Familiar'); end
        text(0.02,0.85,sprintf('n=%d',numel(cellsF)),'Units','normalized','FontSize',8);
        set(gca,'TickDir','out','XTickLabel',[]); box off; xlim([0 1]);
        % --- dF/F top, Novel (col 2) ---
        subplot(nRow,2,rTop+2); hold on;
        cellsN=find(primCat_nov==g);
        drawBand(ctr, tun_nov(cellsN,:), c_nov);
        if gi==1, title('Novel'); end
        text(0.02,0.85,sprintf('n=%d',numel(cellsN)),'Units','normalized','FontSize',8);
        set(gca,'TickDir','out','XTickLabel',[]); box off; xlim([0 1]);
        % --- feature bottom, Familiar (col 1) ---
        subplot(nRow,2,rTop+3); hold on;
        drawFeatureBand(featSig{g}, binIdx, trialIdx, famMask, nBins, ctr, c_feat);
        addRewardLines(rewLocFam);
        ylabel('feature','FontSize',8);
        if gi==nG, xlabel('Track position (norm.)'); end
        set(gca,'TickDir','out'); box off; xlim([0 1]);
        % --- feature bottom, Novel (col 2) ---
        subplot(nRow,2,rTop+4); hold on;
        drawFeatureBand(featSig{g}, binIdx, trialIdx, novMask, nBins, ctr, c_feat);
        addRewardLines(rewLocNov);
        if gi==nG, xlabel('Track position (norm.)'); end
        set(gca,'TickDir','out'); box off; xlim([0 1]);
    end
end

% Mean +/- SEM of a behavioral signal over position bins, computed ACROSS TRIALS.
function drawFeatureBand(sig, binIdx, trialIdx, mask, nBins, ctr, col)
    if isempty(sig), text(0.35,0.5,'(position = axis)','Units','normalized','Color',[.5 .5 .5],'FontSize',8); return; end
    trials = unique(trialIdx(mask & trialIdx>0));
    if isempty(trials), return; end
    perTrial = nan(numel(trials), nBins);
    for ti=1:numel(trials)
        tsel = mask & (trialIdx==trials(ti));
        for b=1:nBins
            s = tsel & (binIdx==b);
            if any(s), perTrial(ti,b)=mean(sig(s),'omitnan'); end
        end
    end
    drawBand(ctr, perTrial, col);   % mean +/- SEM across trials
end

% Distinct normalized reward-zone positions within a context (clustered so
% we draw one line per zone, not one per event).
function locs = rewardPositions(rew, posN, mask)
    locs = [];
    if isempty(rew), return; end
    ev = mask & (rew(:)>0);
    if ~any(ev), return; end
    p = sort(posN(ev));
    % cluster events within 0.03 (normalized) into single zones
    locs = p(1); last = p(1);
    for i=2:numel(p)
        if p(i)-last > 0.03, locs(end+1)=p(i); end %#ok<AGROW>
        last = p(i);
    end
    % collapse near-duplicates to the cluster mean
    if ~isempty(locs)
        keep=true(size(locs));
        for i=2:numel(locs), if locs(i)-locs(i-1)<0.03, keep(i)=false; end, end
        locs=locs(keep);
    end
end

function addRewardLines(locs)
    if isempty(locs), return; end
    yl=ylim;
    for k=1:numel(locs)
        plot([locs(k) locs(k)], yl, '--', 'Color',[0.1 0.6 0.2], 'LineWidth',1);
    end
    ylim(yl);
end

% Reward event vector, mirroring the engine: prefer rewardDs, fall back to
% daqdata.waterValve mapped to imaging frames. Returns a 0/1 (per frame) vector.
function rew = rewardSignal(sData, N)
    rew = zeros(N,1);
    vt = getSig(sData,'rewardDs');
    if ~isempty(vt) && any(double(vt)>0)
        vt=double(vt); rew(1:min(N,numel(vt))) = vt(1:min(N,numel(vt)))>0;
        return;
    end
    % fallback: water valve in daqdata
    if isfield(sData,'daqdata') && isfield(sData.daqdata,'frameIndex') && isfield(sData.daqdata,'waterValve')
        vr=double(sData.daqdata.waterValve>0.5);
        ft=find(diff(sData.daqdata.frameIndex)==1)+1;
        lim=min(N,length(ft)-1);
        for i=1:lim
            if any(vr(ft(i):ft(i+1)-1)), rew(i)=1; end
        end
    end
end

function v = colvec(x, N)
    if isempty(x), v=zeros(N,1); return; end
    x=double(x(:)); if numel(x)>=N, v=x(1:N); else, v=[x; zeros(N-numel(x),1)]; end
end

function s = getSig(sData, name)
    if isfield(sData.behavior.signals, name), s=sData.behavior.signals.(name); else, s=[]; end
end

function p = pupilTraceLocal(sData, N)
    p = nan(N,1);
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

function m = binMean(y, binIdx, mask, nBins)
    m=nan(1,nBins);
    for b=1:nBins
        sel = mask & (binIdx==b);
        if any(sel), m(b)=mean(y(sel),'omitnan'); end
    end
end

function drawBand(x, M, col)
    if isempty(M) || all(isnan(M(:))), text(0.4,0.5,'no cells','Units','normalized','Color',[.5 .5 .5]); return; end
    mu=mean(M,1,'omitnan'); n=sum(~isnan(M),1); sem=std(M,0,1,'omitnan')./sqrt(max(n,1));
    good=~isnan(mu);
    xx=x(good); mu=mu(good); sem=sem(good);
    fill([xx fliplr(xx)],[mu+sem fliplr(mu-sem)],col,'FaceAlpha',0.25,'EdgeColor','none');
    plot(xx,mu,'Color',col,'LineWidth',1.5);
end

function z = zscore_safe(y)
    s=std(y,'omitnan'); if s<1e-9, z=y-mean(y,'omitnan'); else, z=(y-mean(y,'omitnan'))/s; end
end

function s = signCorr(x, y)
    n = min(length(x),length(y));
    r = corr(x(1:n), y(1:n), 'rows','complete');
    if isnan(r), s = 0; else, s = sign(r); if s==0, s=1; end, end
end

% Minimal behavior extraction (velocity + pupil) for the sign test only.
function [vel, pupil_glm, hasPupil, N, famMask, novMask] = quickBehavior(sData, dffLen)
    time_len=length(sData.behavior.signals.timeDs);
    pos_len=length(sData.behavior.signals.corridorPositionDs);
    vel_len=length(sData.behavior.signals.velocityDs);
    N=min([dffLen,time_len,pos_len,vel_len]);
    vel=double(sData.behavior.signals.velocityDs(1:N)); vel=vel(:);
    if length(vel)>2, vel(1:2)=vel(3); end
    pupil=nan(N,1); hasPupil=false;
    try
        if isfield(sData.behavior,'pupildata') && isfield(sData.behavior.pupildata,'diameter')
            pd=double(sData.behavior.pupildata.diameter); pd(isnan(pd))=inf;
            pfi=find(diff(sData.behavior.signals.pupilCamFrameIndexCorr)==1)+1; ifi=find(diff(sData.daqdata.frameIndex)==1);
            pUs=nan(size(sData.behavior.signals.pupilCamFrameIndexCorr)); vl=min(length(pfi),length(pd)); pUs(pfi(1:vl))=pd(1:vl);
            pUs=fillmissing(pUs,'nearest'); pUs(isinf(pUs))=nan; pr=double(pUs(ifi));
            if length(pr)>=N, pupil=pr(1:N); else, pupil=[pr(:);nan(N-length(pr),1)]; end
            pupil=pupil(:); hasPupil=~all(isnan(pupil));
        end
    catch, end
    if hasPupil, pupil_glm=fillmissing(pupil,'nearest'); pupil_glm(isnan(pupil_glm))=mean(pupil_glm,'omitnan');
    else, pupil_glm=zeros(N,1); end

    % context masks (same logic as vipGLM_fitCell) for within-context sign tests
    ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2);
    ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
    trialIdx=sData.behavior.signals.trialIndexesDs(1:N); trialIdx=trialIdx(:);
    vtr=trialIdx>0 & trialIdx<=length(ctx);
    cctx=nan(N,1); cctx(vtr)=ctx(trialIdx(vtr)); cctx=fillmissing(cctx,'nearest');
    fam_id=ctx(1); uctx=unique(ctx); nov_id=uctx(uctx~=fam_id);
    famMask=(cctx==fam_id);   % all familiar frames (both blocks)
    if isempty(nov_id), novMask=false(N,1); else, novMask=(cctx==nov_id(1)); end
end

function closeIfValid(d), if isvalid(d), close(d); end, end

function progressTick(d, n, doReset)
    persistent cnt
    if isempty(cnt)||doReset, cnt=0; end
    if doReset, return; end
    cnt=cnt+1;
    if isvalid(d), d.Value=0.02+0.93*(cnt/n); d.Message=sprintf('Fitting cell %d/%d...',cnt,n); end
end
