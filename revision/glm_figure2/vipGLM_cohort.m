function vipGLM_cohort()
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData');
    if ~exist(dfSess,'dir'), dfSess=userpath; end
% VIPGLM_COHORT  Cohort-level population GLM across multiple sessions -- the main entry point.
%
% AIM: run the per-cell encoding model (vipGLM_fitCell) on every cell of every
% selected session and aggregate into one results table + summary figures.
%
% INPUT:  the SESSION DATA files to run -- ONE OR MORE sData .mat (multi-select prompt).
% OUTPUT: saves cohortGLM_<n>sessions_<date>.mat (res.cellTable + stats) and shows
%         summary figures. Run vipGLM_paperFigure afterwards for the Fig.1 strip.
% NEEDS:  vipGLM_fitCell.m on the path.
%
%   Aggregates PER SESSION (session = unit of observation); counts are expressed
%   as FRACTIONS of cells so sessions with different cell numbers are comparable.
%
%   No inferential statistics are computed here: the cohort saves the per-cell
%   cellTable + aggregate means/SEM and draws the summary figures only. Any
%   significance testing is done in dedicated standalone scripts (e.g.
%   vipGLM_compareCustomGroups, vipGLM_clusterValidation) run on the saved .mat.
%
%   Requires vipGLM_fitCell.m on the path.

    uiwait(msgbox({'Select the SESSION DATA files to run.','',...
        'ONE OR MORE m8xxx session .mat files (hold Cmd/Shift for several).'},'Session data (one or more)','help'));
    [files, pathName] = uigetfile('*.mat','SESSION DATA files (ONE OR MORE m8xxx .mat; multi-select)', ...
        dfSess, 'MultiSelect','on');
    if isequal(files,0), return; end
    if ischar(files), files={files}; end           % single selection -> cell
    nSess = numel(files);

    uiFig = uifigure('Name','Cohort GLM','Position',[100 100 380 490]);
    g = uigridlayout(uiFig,[12 2]); g.RowHeight=repmat({'fit'},1,12); g.ColumnWidth={170,'1x'};
    uilabel(g,'Text','Ridge \lambda:');        lamF  = uieditfield(g,'numeric','Value',0.5);
    uilabel(g,'Text','Lambda mode:');           lamDD = uidropdown(g,'Items',{'Fixed lambda','CV-tuned lambda'},'Value','Fixed lambda');
    uilabel(g,'Text','Significance:');          sigDD = uidropdown(g,'Items',{'Permutation (rigorous)','t-test (fast)'},'Value','Permutation (rigorous)');
    uilabel(g,'Text','# permutations:');        permF = uieditfield(g,'numeric','Value',200,'Limits',[20 2000],'RoundFractionalValues','on');
    uilabel(g,'Text','Lag min (s):');           lagMinF = uieditfield(g,'numeric','Value',-3,'Limits',[-10 10]);
    uilabel(g,'Text','Lag max (s):');           lagMaxF = uieditfield(g,'numeric','Value',3,'Limits',[-10 10]);
    uilabel(g,'Text','# lag points:');          nLagF   = uieditfield(g,'numeric','Value',5,'Limits',[1 21],'RoundFractionalValues','on');
    uilabel(g,'Text','Position basis:');        posTypeDD = uidropdown(g,'Items',{'Binary bins (Turi)','Gaussian RBF'},'Value','Binary bins (Turi)');
    uilabel(g,'Text','# position bins (0=off):'); nPosF = uieditfield(g,'numeric','Value',100,'Limits',[0 200],'RoundFractionalValues','on');
    uilabel(g,'Text','Use parfor (cells):');    parDD = uidropdown(g,'Items',{'Yes','No'},'Value','Yes');
    uilabel(g,'Text','Save results .mat:');     saveDD= uidropdown(g,'Items',{'Yes','No'},'Value','Yes');
    runB = uibutton(g,'Text',sprintf('Run cohort (%d sessions)',nSess),'ButtonPushedFcn', @(s,e) runCohort( ...
        files, pathName, lamF.Value, round(permF.Value), strcmp(parDD.Value,'Yes'), strcmp(saveDD.Value,'Yes'), ...
        ternary(contains(lamDD.Value,'CV'),'nested','fixed'), ternary(contains(sigDD.Value,'t-test'),'ttest','perm'), ...
        lagMinF.Value, lagMaxF.Value, round(nLagF.Value), round(nPosF.Value), ...
        ternary(contains(posTypeDD.Value,'RBF'),'rbf','binary'), uiFig));
    runB.Layout.Row=12; runB.Layout.Column=[1 2];
end

function v = ternary(cond,a,b), if cond, v=a; else, v=b; end, end

% One row of the per-cell results table (session x cell x context). Features in
% engine order: 1 Vel 2 Accel 3 Lick 4 Reward 5 Position 6 Pupil.
function r = makeCellRow(sessName, sessIdx, cellID, o, context, actFam, actNov, velSign)
    if nargin<8 || isempty(velSign), velSign=1; end
    fn={'Velocity','Acceleration','Licking','Reward','Position','Pupil'};
    r.session=string(sessName); r.sessIdx=sessIdx; r.cellID=cellID; r.context=string(context);
    r.actMean_fam=actFam; r.actMean_nov=actNov; r.actDelta=actNov-actFam;  % activity change novelty-familiar
    % +1/-1 velocity modulation sign; NaN preserved (undefined: dead/flat cell)
    if isnan(velSign), r.velSign=NaN; else, r.velSign=sign(velSign); if r.velSign==0, r.velSign=1; end, end
    if isempty(o) || ~isfield(o,'ok') || ~o.ok
        r.ok=false; r.R2_full=NaN; r.hasPupil=false; r.pupilDetFrac=NaN; r.lambdaUsed=NaN; r.nFramesUsed=NaN;
        for f=1:6
            r.(['dR2_' fn{f}])=NaN; r.(['p_' fn{f}])=NaN; r.(['sig_' fn{f}])=false; r.(['singleR2_' fn{f}])=NaN;
        end
        r.nSigFeatures=NaN; r.primaryFeature="none"; r.primaryFeatureSigned="none";
        return;
    end
    r.ok=true; r.R2_full=o.R2_full; r.hasPupil=o.hasPupil;
    r.pupilDetFrac=getfielddef(o,'pupilDetFrac',NaN); r.lambdaUsed=o.lambdaUsed; r.nFramesUsed=o.nFramesUsed;
    for f=1:6
        r.(['dR2_' fn{f}])      = o.dR2_unique(f);
        r.(['p_' fn{f}])        = o.p_vals(f);
        r.(['sig_' fn{f}])      = logical(o.sig(f));
        r.(['singleR2_' fn{f}]) = o.single_R2(f);
    end
    r.nSigFeatures = sum(o.sig & ~isnan(o.dR2_unique));
    valid=find(o.sig & ~isnan(o.dR2_unique));
    if isempty(valid)
        r.primaryFeature="none"; r.primaryFeatureSigned="none";
    else
        [~,m]=max(o.dR2_unique(valid)); pf=fn{valid(m)};
        r.primaryFeature=string(pf);
        if strcmp(pf,'Velocity')   % signed label: Velocity(+) / Velocity(-) / Velocity(?) if undefined
            if isnan(r.velSign), sgn='?'; else, sgn=tern(r.velSign>=0,'+','-'); end
            r.primaryFeatureSigned=string(sprintf('Velocity(%s)', sgn));
        else
            r.primaryFeatureSigned=string(pf);
        end
    end
end

function s = tern(c,a,b), if c, s=a; else, s=b; end, end

% Save every open figure as a PNG into a timestamped folder on the Desktop.
function stampAllFigures(paramStr)
% Add a small parameter footnote to every open figure (task 2), so saved PNGs
% are self-documenting for parametric exploration.
    figs = findobj('Type','figure');
    for k=1:numel(figs)
        f=figs(k);
        try
            annotation(f,'textbox',[0 0 1 0.03],'String',paramStr, ...
                'FontSize',7,'Color',[0.35 0.35 0.35],'EdgeColor','none', ...
                'HorizontalAlignment','center','VerticalAlignment','bottom', ...
                'Interpreter','tex','FitBoxToText','off');
        catch
        end
    end
end

function saveFiguresToDesktop()
    % locate Desktop cross-platform, with fallback to the home dir / pwd
    home = char(java.lang.System.getProperty('user.home'));
    if isempty(home), home = pwd; end
    desktop = fullfile(home,'Desktop');
    if ~isfolder(desktop), desktop = home; end          % fall back if no Desktop
    folder = fullfile(desktop, ['vipGLM_figures_' datestr(now,'yyyy-mm-dd_HHMM')]);
    [ok,~] = mkdir(folder);
    if ~ok
        warning('Could not create figure folder at %s; figures not saved.', folder); return;
    end
    figs = findobj('Type','figure');
    figs = flipud(figs(:));                              % save in creation order
    n = 0;
    for i = 1:numel(figs)
        f = figs(i);
        nm = get(f,'Name');
        if isempty(nm), nm = sprintf('figure_%d', double(f.Number)); end
        nm = regexprep(nm, '[^\w\-]+', '_');             % filesystem-safe
        nm = regexprep(nm, '_+', '_'); nm = regexprep(nm, '^_|_$', '');
        fname = fullfile(folder, sprintf('%02d_%s.png', i, nm));
        try
            exportgraphics(f, fname, 'Resolution', 200); % crisp PNG (R2020a+)
            n = n + 1;
        catch
            try, saveas(f, fname); n = n + 1;            % fallback for older MATLAB
            catch, warning('Could not save figure %d', i); end
        end
    end
    fprintf('Saved %d figure(s) to: %s\n', n, folder);
end

function r = emptyCellRow()
    r = makeCellRow("",0,0,[],"",NaN,NaN); r(1)=[];   % 0x1 struct with the right fields
end

function v = getfielddef(s,f,def), if isfield(s,f), v=s.(f); else, v=def; end, end

% Keep only the fields vipGLM_fitCell actually reads, so parfor broadcasts a
% tiny struct instead of the full session (raw imaging). imdata is NOT needed
% by the engine (it indexes the dff_all argument directly).
function sd = trimSData(sData)
    sd = struct();
    % behavior.signals (the vectors the engine extracts)
    bs = sData.behavior.signals;
    keepSig = {'timeDs','corridorPositionDs','velocityDs','lickIfreqDs', ...
               'rewardDs','trialIndexesDs','pupilCamFrameIndexCorr'};
    s2 = struct();
    for k=1:numel(keepSig)
        % case-insensitive match: copy whatever case the source uses into the
        % canonical (lowercase) name the engine reads. Guards against e.g.
        % 'VelocityDs' in some pipelines vs 'velocityDs' expected here.
        src = matchFieldCI(bs, keepSig{k});
        if ~isempty(src), s2.(keepSig{k}) = bs.(src); end
    end
    sd.behavior.signals = s2;
    % trial matrices (contextMatrix)
    if isfield(sData.behavior,'trialMatrices') && isfield(sData.behavior.trialMatrices,'contextMatrix')
        sd.behavior.trialMatrices.contextMatrix = sData.behavior.trialMatrices.contextMatrix;
    end
    % pupil
    if isfield(sData.behavior,'pupildata')
        sd.behavior.pupildata = sData.behavior.pupildata;
    end
    % daqdata (frameIndex + waterValve, used in the reward/pupil fallback paths)
    if isfield(sData,'daqdata')
        dd = struct();
        if isfield(sData.daqdata,'frameIndex'), dd.frameIndex=sData.daqdata.frameIndex; end
        if isfield(sData.daqdata,'waterValve'), dd.waterValve=sData.daqdata.waterValve; end
        sd.daqdata = dd;
    end
    % NOTE: imdata is intentionally omitted -- vipGLM_fitCell reads dF/F from the
    % dff_all argument, never from sData.imdata, so the big imaging struct is not
    % broadcast to the workers.
end

% Case-insensitive fieldname match: returns the actual field name in struct s
% that matches `name` ignoring case, or '' if none.
function fn = matchFieldCI(s, name)
    fn=''; flds=fieldnames(s);
    idx=find(strcmpi(flds,name),1);
    if ~isempty(idx), fn=flds{idx}; end
end

% =====================================================================
function runCohort(files, pathName, lambdaVal, nPerms, usePar, doSave, lambdaMode, sigTest, lagMin, lagMax, nLags, nPosRBF, posBasis, uiFig)
    tStart = tic;   % wall-clock timing for the whole cohort run
    nSess=numel(files);
    nCat=7; nMux=7;   % 7 categories (velocity sign-split, pupil single); multiplexing bins 0..6
    cat_names={'Vel (+)','Vel (-)','Accel','Licking','Reward','Position','Pupil'};

    % per-session aggregates (fractions). rows = sessions, cols = category/bin.
    sigFrac_fam=nan(nSess,nCat); sigFrac_nov=nan(nSess,nCat);
    primFrac_fam=nan(nSess,nCat); primFrac_nov=nan(nSess,nCat);
    muxFrac_fam=nan(nSess,nMux); muxFrac_nov=nan(nSess,nMux);
    dR2_fam=nan(nSess,nCat); dR2_nov=nan(nSess,nCat);   % mean dR2 over significant cells
    nCellsPer=nan(nSess,1); sessNames=cell(nSess,1);
    pupilDet=nan(nSess,1);            % per-session pupil detection fraction (diagnostic)
    sessTime=nan(nSess,1);            % per-session wall-clock time (s)
    cellRows=emptyCellRow();          % growing struct array, one row per cell x context
    % cross-session pools for the grouped FEATURE-tuning figure
    GT.nBins=12; GT.nG=7;   % feature binned into quantiles -> fewer bins than position
    GT.tunFam=cell(1,GT.nG); GT.tunNov=cell(1,GT.nG);   % pooled per-cell dF/F vs feature
    GT.xcFam =cell(1,GT.nG); GT.xcNov =cell(1,GT.nG);   % matching feature-value bin centers
    GT.featName={'Velocity','Velocity','Acceleration','Lick rate','Reward','Position','Pupil diameter'};

    % cross-session pools for the STACKED position-tuning figure (dF/F vs track
    % position + defining feature underneath + reward lines), grouped by winning feature.
    PT.nBins=40; PT.nG=7;
    PT.tunFam=cell(1,PT.nG); PT.tunNov=cell(1,PT.nG);    % pooled per-cell dF/F vs position
    PT.featFam=cell(1,PT.nG); PT.featNov=cell(1,PT.nG);  % per-session feature curves vs position
    PT.rewFam={}; PT.rewNov={};                           % per-session reward locations

    d = uiprogressdlg(uiFig,'Title','Cohort GLM','Message','Starting...','Value',0);
    cleanup = onCleanup(@() closeIfValid(d)); %#ok<NASGU>

    for si=1:nSess
        if ~isvalid(d), fprintf('User canceled the cohort run.\n'); return; end
        fn=files{si}; sessNames{si}=fn;
        d.Value=(si-1)/nSess; d.Message=sprintf('Session %d/%d: %s',si,nSess,fn);
        fprintf('Session %d/%d: %s\n', si, nSess, fn);   % console trail (independent of GUI dialog)
        tSess=tic;
        try
            S=load(fullfile(pathName,fn)); sData=S.sData;
        catch ME
            warning('Skipping %s: %s',fn,ME.message); continue;
        end
        if ~(isfield(sData.imdata,'roiSignals') && numel(sData.imdata.roiSignals)>=2)
            warning('Skipping %s: no VIP roiSignals(2).dff',fn); continue;
        end
        dff_all=sData.imdata.roiSignals(2).dff;
        nCells=size(dff_all,1); nCellsPer(si)=nCells;

        % within-context sign masks (same logic as the engine)
        [velTrace,pupilTrace,~,N,famMask,novMask]=quickBehavior(sData,size(dff_all,2));

        % Build a TRIMMED sData with only the fields the engine reads. The full
        % struct (raw imaging) is huge; broadcasting it to parfor workers every
        % session is what makes later sessions progressively slower. The trimmed
        % copy is tiny (a few behavior vectors), so the per-session broadcast cost
        % stays flat. dff_all is passed separately and indexed directly.
        sd = trimSData(sData);

        % fit all cells, both contexts
        famRes=cell(nCells,1); novRes=cell(nCells,1); seeds=1:nCells;
        runP = usePar && ~isempty(ver('parallel'));
        if runP
            dq = parallel.pool.DataQueue;
            progressTick(d, nCells, true);                       % reset cell counter
            afterEach(dq, @(~) progressTick(d, nCells, false));   % tick on each worker report
            parfor c=1:nCells
                rng(seeds(c));
                dffc = dff_all(c,:);   % sliced: only this row is sent to the worker
                famRes{c}=vipGLM_fitCell(sd,dffc,1,'Familiar Context',lambdaVal,nPerms,lambdaMode,sigTest,lagMin,lagMax,nLags,nPosRBF,posBasis);
                novRes{c}=vipGLM_fitCell(sd,dffc,1,'Novel Context',   lambdaVal,nPerms,lambdaMode,sigTest,lagMin,lagMax,nLags,nPosRBF,posBasis);
                send(dq, c);
            end
        else
            for c=1:nCells
                rng(seeds(c));
                dffc = dff_all(c,:);   % slice (match parfor branch)
                famRes{c}=vipGLM_fitCell(sd,dffc,1,'Familiar Context',lambdaVal,nPerms,lambdaMode,sigTest,lagMin,lagMax,nLags,nPosRBF,posBasis);
                novRes{c}=vipGLM_fitCell(sd,dffc,1,'Novel Context',   lambdaVal,nPerms,lambdaMode,sigTest,lagMin,lagMax,nLags,nPosRBF,posBasis);
                if isvalid(d)
                    d.Value=(si-1)/nSess + (c/nCells)*(1/nSess);
                    d.Message=sprintf('Session %d/%d: cell %d/%d', si, nSess, c, nCells);
                else
                    fprintf('User canceled the cohort run.\n'); return;
                end
            end
        end

        % aggregate this session
        [sigFrac_fam(si,:),primFrac_fam(si,:),muxFrac_fam(si,:),dR2_fam(si,:)] = ...
            aggregateSession(famRes, dff_all, velTrace, pupilTrace, famMask, N);
        [sigFrac_nov(si,:),primFrac_nov(si,:),muxFrac_nov(si,:),dR2_nov(si,:)] = ...
            aggregateSession(novRes, dff_all, velTrace, pupilTrace, novMask, N);

        % ---- per-cell record (one row per cell x context) for later re-analysis ----
        % captures every feature's dR2/p/significance, plus simple activity measures
        % per context so questions like "pupil-tuned cells that gained activity in
        % novelty" can be answered from the saved table without re-running.
        % per-cell activity means + velocity sign, written into the saved cellTable
        for c=1:nCells
            yfull = double(dff_all(c,1:N))';
            actFam = mean(yfull(famMask),'omitnan');
            actNov = mean(yfull(novMask),'omitnan');
            vSignF = signCorr(velTrace(famMask), yfull(famMask));   % velocity sign within familiar
            vSignN = signCorr(velTrace(novMask), yfull(novMask));   % velocity sign within novel
            cellRows(end+1) = makeCellRow(fn, si, c, famRes{c}, 'Familiar', actFam, actNov, vSignF); %#ok<AGROW>
            cellRows(end+1) = makeCellRow(fn, si, c, novRes{c}, 'Novel',    actFam, actNov, vSignN); %#ok<AGROW>
        end
        if ~isempty(famRes) && ~isempty(famRes{1}) && isfield(famRes{1},'pupilDetFrac')
            pupilDet(si) = famRes{1}.pupilDetFrac;
        end

        % accumulate grouped position-tuning data (pooled cells; per-session features)
        GT = accumGroupTuning(GT, dff_all, sData, famRes, novRes, velTrace, pupilTrace, famMask, novMask, N);

        % accumulate stacked position-tuning data (dF/F vs position + feature underlay)
        PT = accumPositionTuning(PT, dff_all, sData, famRes, velTrace, pupilTrace, famMask, novMask, N);

        % free large per-session variables so memory does not accumulate
        clear S sData sd dff_all famRes novRes velTrace pupilTrace

        secSess=toc(tSess); sessTime(si)=secSess;
        fprintf('   session %d done in %.1f s (%d cells; cumulative %.1f s)\n', ...
            si, secSess, nCellsPer(si), sum(sessTime,'omitnan'));
    end
    d.Value=0.95; d.Message='Aggregating + statistics...';

    % drop sessions that failed entirely (all NaN row in sig fraction)
    good = ~all(isnan(sigFrac_fam),2);
    sigFrac_fam=sigFrac_fam(good,:); sigFrac_nov=sigFrac_nov(good,:);
    primFrac_fam=primFrac_fam(good,:); primFrac_nov=primFrac_nov(good,:);
    muxFrac_fam=muxFrac_fam(good,:); muxFrac_nov=muxFrac_nov(good,:);
    dR2_fam=dR2_fam(good,:); dR2_nov=dR2_nov(good,:);
    nGood=sum(good);

    % ---- figures: 4-panel population summary (no omnibus stats annotations) ----
    plotCohort(sigFrac_fam,sigFrac_nov,primFrac_fam,primFrac_nov, ...
        muxFrac_fam,muxFrac_nov,dR2_fam,dR2_nov,cat_names,nGood);

    % grouped feature-tuning figure (pooled across all cells of all sessions)
    plotCohortGroupedTuning(GT, cat_names);

    % stacked position tuning: dF/F vs position + defining feature + reward lines
    plotCohortPositionTuning(PT, cat_names);

    % param string stamped on every figure
    paramStr = sprintf(['\\lambda=%.3g (%s) | sig=%s%s | lags %.2g:%.2g s (%d pts), \\sigma=0.125s | ' ...
        'pos: %d %s bins | parfor=%d | N=%d sessions'], ...
        lambdaVal, lambdaMode, sigTest, ternStr(strcmpi(sigTest,'perm'),sprintf(', perms=%d',nPerms),''), ...
        lagMin, lagMax, nLags, nPosRBF, posBasis, usePar, nGood);

    if doSave
        outFile=sprintf('cohortGLM_%dsessions_%s.mat',nGood,datestr(now,'yyyymmdd_HHMMSS'));
        % flat per-cell table for later re-analysis (filter/join without re-running).
        % This table + meta is what every downstream script (paperFigure, etc.) reads.
        if ~isempty(cellRows), cellTable=struct2table(cellRows); else, cellTable=table(); end
        res=struct( ...
            'sigFrac_fam',sigFrac_fam,'sigFrac_nov',sigFrac_nov, ...
            'primFrac_fam',primFrac_fam,'primFrac_nov',primFrac_nov, ...
            'muxFrac_fam',muxFrac_fam,'muxFrac_nov',muxFrac_nov, ...
            'dR2_fam',dR2_fam,'dR2_nov',dR2_nov,'cat_names',{cat_names}, ...
            'cellTable',{cellTable},'pupilDetFrac',pupilDet(good), ...
            'nSessions',nGood,'nCellsPer',nCellsPer(good),'sessNames',{sessNames(good)}, ...
            'meta',struct( ...
                'lambda',lambdaVal,'lambdaMode',lambdaMode,'sigTest',sigTest,'nPerms',nPerms, ...
                'useParfor',usePar,'sessionPath',pathName,'sessionFiles',{{files}}, ...
                'familiarDefinition','all familiar frames (both blocks)', ...
                'shifts_s',linspace(lagMin,lagMax,nLags),'lagMin_s',lagMin,'lagMax_s',lagMax,'nLags',nLags,'kernel_sigma_s',0.125,'nPosRBF',nPosRBF,'posBasis',posBasis, ...
                'cv','blocked 10s chunks, 10-fold','rngSeedPerCell','1:nCells', ...
                'pupilGuard','none (pupilDetFrac recorded as diagnostic only)', ...
                'engine','vipGLM_fitCell','runDate',datestr(now,'yyyy-mm-dd HH:MM:SS'), ...
                'matlabVersion',version)); %#ok<NASGU>
        save(outFile,'res'); fprintf('Saved %s\n',outFile);
    end

    % console summary
    fprintf('\n=== Cohort: %d sessions ===\n',nGood);
    pd=pupilDet(good);
    fprintf('Pupil detection fraction per session (threshold 0.5 -> pupil used):\n');
    fprintf('  %s\n', num2str(pd','%.2f  '));
    fprintf('  sessions with usable pupil: %d/%d\n', sum(pd>=0.5 & ~isnan(pd)), numel(pd));
    close(d);

    % ---- stamp every open figure with the parameter string (task 2) ----
    stampAllFigures(paramStr);

    % ---- save all open figures as PNG to a timestamped Desktop folder ----
    saveFiguresToDesktop();

    elapsed = toc(tStart);
    h=floor(elapsed/3600); m=floor(mod(elapsed,3600)/60); s=mod(elapsed,60);
    fprintf('\n=== Cohort run completed in %dh %dm %.1fs (%.1f s total) ===\n', h, m, s, elapsed);
    if any(~isnan(sessTime))
        fprintf('    per-session: mean %.1f s, min %.1f s, max %.1f s (n=%d sessions timed)\n', ...
            mean(sessTime,'omitnan'), min(sessTime), max(sessTime), sum(~isnan(sessTime)));
    end
end

function [sigFrac,primFrac,muxFrac,dR2mean,dR2mean6] = aggregateSession(res, dff_all, velTrace, pupilTrace, mask, N)
    nCat=7; nMux=7;
    sigCount=zeros(1,nCat); primCount=zeros(1,nCat); dR2acc=cell(1,nCat);
    muxCount=zeros(1,nMux); nUsed=0;
    dR2acc6=cell(1,6);   % sign-agnostic: 1 Vel 2 Accel 3 Lick 4 Reward 5 Position 6 Pupil
    map=[3 4 5 6];   % accel,lick,reward,position -> categories 3..6
    sessHasPupil=false;   % becomes true if any fitted cell had a working pupil camera

    for c=1:numel(res)
        o=res{c};
        if isempty(o)||~o.ok, continue; end
        nUsed=nUsed+1;
        if o.hasPupil, sessHasPupil=true; end
        sig=o.sig; dR2=o.dR2_unique;

        % sign-agnostic 6-category dR2 (velocity & pupil NOT split by sign)
        for f=1:6
            if sig(f) && ~isnan(dR2(f)), dR2acc6{f}=[dR2acc6{f}, dR2(f)]; end
        end

        % multiplexing: # significant features (pupil NaN-safe)
        nFeat=sum(sig & ~isnan(dR2));
        muxCount(min(nFeat,6)+1)=muxCount(min(nFeat,6)+1)+1;

        % primary = max unique dR2 among significant
        primIdx=0; valid=find(sig & ~isnan(dR2));
        if ~isempty(valid), [~,m]=max(dR2(valid)); primIdx=valid(m); end

        % within-context sign for vel/pupil
        ydff=double(dff_all(c,1:N))';
        vS=signCorr(velTrace(mask), ydff(mask));
        pS=signCorr(pupilTrace(mask), ydff(mask));

        % velocity -> cat 1(+)/2(-); skip if sign undefined (dead/flat cell)
        if sig(1) && ~isnan(vS)
            idx=1*(vS>=0)+2*(vS<0);
            sigCount(idx)=sigCount(idx)+1; dR2acc{idx}=[dR2acc{idx},dR2(1)];
            if primIdx==1, primCount(idx)=primCount(idx)+1; end
        end
        % accel/lick/reward/position -> 3..6
        for f=2:5
            if sig(f)
                idx=map(f-1);
                sigCount(idx)=sigCount(idx)+1; dR2acc{idx}=[dR2acc{idx},dR2(f)];
                if primIdx==f, primCount(idx)=primCount(idx)+1; end
            end
        end
        % pupil -> 7 (single, no sign split)
        if sig(6) && ~isnan(dR2(6))
            idx=7;
            sigCount(idx)=sigCount(idx)+1; dR2acc{idx}=[dR2acc{idx},dR2(6)];
            if primIdx==6, primCount(idx)=primCount(idx)+1; end
        end
    end

    if nUsed==0
        % UNMASK SILENT FAILURES: report WHY every cell was rejected. Tally the
        % engine's msg across all cells (cell 1 may differ from the rest), so a
        % whole-session drop (e.g. missing reward vector) is visible, not silent.
        reasons = strings(0,1);
        for q=1:numel(res)
            o=res{q};
            if ~isempty(o) && isfield(o,'msg') && ~isempty(o.msg), reasons(end+1,1)=string(o.msg); end %#ok<AGROW>
        end
        if isempty(reasons)
            fprintf('      [!] Session dropped (all %d cells rejected). Reason: unknown\n', numel(res));
        else
            u=unique(reasons); cnt=arrayfun(@(r) sum(reasons==r), u);
            [~,o]=sort(cnt,'descend');
            msg=strjoin(arrayfun(@(i) sprintf('%s (%dx)',u(o(i)),cnt(o(i))), 1:numel(u), 'uni',0), ', ');
            fprintf('      [!] Session dropped (%d/%d cells rejected). Reason(s): %s\n', sum(cnt), numel(res), msg);
        end
        sigFrac=nan(1,nCat); primFrac=nan(1,nCat); muxFrac=nan(1,nMux); dR2mean=nan(1,nCat); dR2mean6=nan(1,6); return;
    end
    sigFrac = sigCount / nUsed;     % fraction of cells significant per category
    primFrac= primCount / nUsed;    % fraction whose primary feature is this category
    muxFrac = muxCount / nUsed;     % fraction of cells with k significant features
    dR2mean = cellfun(@(v) mean(v,'omitnan'), dR2acc);   % mean dR2 among SIGNIFICANT cells
    % Sessions with no significant cells for a category stay NaN (not 0): this is
    % a tuning-STRENGTH measure conditioned on significance, so "no tuned cells"
    % means "no strength data", not "zero strength". Cross-session averaging then
    % uses 'omitnan' so empty sessions don't dilute the strength estimate.
    dR2mean(cellfun(@isempty,dR2acc)) = NaN;
    dR2mean6 = cellfun(@(v) mean(v,'omitnan'), dR2acc6);
    dR2mean6(cellfun(@isempty,dR2acc6)) = NaN;   % no sig cells -> no strength data (not 0)

    % MISSING-HARDWARE GUARD: if this session had no pupil camera, the pupil
    % categories are "no data", not "zero effect". Set them to NaN so they are
    % dropped from cross-session means/SEM and from the repeated-measures stats
    % (which use complete cases) rather than dragging pupil toward zero.
    if ~sessHasPupil
        sigFrac(7)=NaN; primFrac(7)=NaN; dR2mean(7)=NaN;   % pupil (single category)
        dR2mean6(6)=NaN;                                               % merged pupil
    end
end

function padj = holm(p)
    padj=nan(size(p)); valid=~isnan(p); pv=p(valid);
    [s,ord]=sort(pv); m=numel(pv); adj=zeros(1,m); run=0;
    for i=1:m
        run=max(run,(m-i+1)*s(i)); adj(i)=min(run,1);
    end
    out=nan(1,m); out(ord)=adj; padj(valid)=out;
end

% =====================================================================
function plotCohort(sigF,sigN,primF,primN,muxF,muxN,dR2F,dR2N,cat_names,nS)
    f=figure('Color','w','Position',[80 80 1500 850]); %#ok<NASGU>
    sgtitle(sprintf('Cohort population encoding  (mean \\pm SEM, N=%d sessions)',nS),'FontWeight','bold','FontSize',14);
    c_fam=[0.2 0.5 0.8]; c_nov=[0.8 0.3 0.3]; nCat=numel(cat_names);

    msem=@(M) deal(mean(M,1,'omitnan'), std(M,0,1,'omitnan')./sqrt(sum(~isnan(M),1)));

    % Panel 1: significant fraction
    subplot(2,2,1); hold on;
    [m1f,e1f]=msem(sigF); [m1n,e1n]=msem(sigN);
    groupedBarSEM(1:nCat,[m1f;m1n]',[e1f;e1n]',c_fam,c_nov);
    set(gca,'XTick',1:nCat,'XTickLabel',cat_names,'XTickLabelRotation',30,'TickDir','out');
    ylabel('Fraction of cells'); title('1. Significant modulation (fraction)'); box off;
    legend({'Familiar','Novel'},'Location','best');

    % Panel 2: primary feature fraction
    subplot(2,2,2); hold on;
    [m2f,e2f]=msem(primF); [m2n,e2n]=msem(primN);
    groupedBarSEM(1:nCat,[m2f;m2n]',[e2f;e2n]',c_fam,c_nov);
    set(gca,'XTick',1:nCat,'XTickLabel',cat_names,'XTickLabelRotation',30,'TickDir','out');
    ylabel('Fraction of cells'); title('2. Primary (winning) feature (fraction)'); box off;

    % Panel 3: multiplexing distribution
    subplot(2,2,3); hold on;
    [m3f,e3f]=msem(muxF); [m3n,e3n]=msem(muxN);
    groupedBarSEM(0:6,[m3f;m3n]',[e3f;e3n]',c_fam,c_nov);
    set(gca,'XTick',0:6,'TickDir','out');
    xlabel('# significant features per cell'); ylabel('Fraction of cells');
    title('3. Multiplexing distribution'); box off;

    % Panel 4: effect size dR2
    subplot(2,2,4); hold on;
    [m4f,e4f]=msem(dR2F); [m4n,e4n]=msem(dR2N);
    groupedBarSEM(1:nCat,[m4f;m4n]',[e4f;e4n]',c_fam,c_nov);
    set(gca,'XTick',1:nCat,'XTickLabel',cat_names,'XTickLabelRotation',30,'TickDir','out');
    ylabel('Mean unique \DeltaR^2'); title('4. Effect size'); box off;
end

function drawOnePmat(pos, P, cat_names, ttl)
    nC=numel(cat_names);
    subplot(1,2,pos);
    L = -log10(P);                       % larger = smaller p
    Lmax = max([L(~isnan(L)); 3]);       % cap reference for contrast/scale
    Ld = L; Ld(isnan(Ld))=0;
    imagesc(Ld,[0 Lmax]); axis square; colormap(gca, parula);  % parula: dark=low, bright=high
    set(gca,'XTick',1:nC,'XTickLabel',cat_names,'XTickLabelRotation',30, ...
        'YTick',1:nC,'YTickLabel',cat_names,'TickDir','out');
    cb=colorbar; cb.Label.String='-log_{10}(p)';
    % mark the p=0.05 line on the colorbar (-log10(0.05)=1.30)
    thr=-log10(0.05);
    if thr<=Lmax
        cb.Ticks=sort(unique([0 thr round(Lmax)]));
        lbl=arrayfun(@(t) sprintf('%.2g',t),cb.Ticks,'uni',0);
        [~,it]=min(abs(cb.Ticks-thr)); lbl{it}=sprintf('%.2f (p=.05)',thr);
        cb.TickLabels=lbl;
    end
    title(ttl);
    hold on;
    for a=1:nC
        for b=1:nC
            if a==b, continue; end
            p=P(a,b); if isnan(p), continue; end
            if p<0.05
                rectangle('Position',[b-0.5 a-0.5 1 1],'EdgeColor',[1 0 0],'LineWidth',2);
            end
            txt=sprintf('%.3f',p); if p<0.001, txt='<.001'; end
            % text color by cell brightness: parula is dark at low values (white text)
            % and bright/yellow at high values (black text) -> always readable.
            frac=Ld(a,b)/max(Lmax,eps);
            lc=[1 1 1]; if frac>0.55, lc=[0 0 0]; end
            text(b,a,txt,'HorizontalAlignment','center','FontSize',7,'Color',lc,'FontWeight','bold');
        end
    end
    hold off;
end


function s = sigStars(p)
    if p<0.001, s='***'; elseif p<0.01, s='**'; elseif p<0.05, s='*'; else, s='n.s.'; end
end

function groupedBarSEM(xc, M, E, c1, c2)
    b=bar(xc, M, 'grouped', 'EdgeColor','k');
    b(1).FaceColor=c1; b(2).FaceColor=c2;
    for i=1:2
        errorbar(b(i).XEndPoints, b(i).YData, E(:,i)', 'k', 'LineStyle','none','LineWidth',1);
    end
end

% =====================================================================
function s = signCorr(x, y)
    n=min(length(x),length(y));
    r=corr(x(1:n),y(1:n),'rows','complete');
    % NaN (flat dF/F, no behaviour variance, dead cell) -> NaN, NOT a default
    % sign. Callers treat NaN as "no sign" and skip the +/- split, so inactive
    % cells don't inflate the (+) categories.
    if isnan(r), s=NaN; else, s=sign(r); if s==0, s=1; end, end
end

% =====================================================================
% Accumulate grouped FEATURE-tuning across sessions: for each group, bin each
% cell's dF/F by the VALUE of that group's defining feature (velocity for the
% velocity groups, pupil diameter for pupil, etc.), pooled across all cells of
% all sessions. This shows what each group actually encodes -- unlike binning
% everything by position, which just recovers the shared position structure.
% Feature is binned into nBins quantile bins (robust to the very different
% scales/distributions of velocity vs pupil vs position).
function GT = accumGroupTuning(GT, dff_all, sData, famRes, novRes, velTrace, pupilTrace, famMask, novMask, N)
    nBins=GT.nBins;
    pos=double(sData.behavior.signals.corridorPositionDs(1:N)); pos=pos(:);
    vel=velTrace(:);
    timeA=double(sData.behavior.signals.timeDs(1:N)); FsA=1/mean(diff(timeA),'omitnan');
    acc=gradient(vel)*FsA;   % match GLM engine (cm/s^2), not [0;diff(vel)]
    lick=colvecC(getSigC(sData,'lickIfreqDs'),N); pup=pupilTrace(:);
    % group -> defining feature signal (1 Vel+ 2 Vel- 3 Accel 4 Lick 5 Reward 6 Position 7 Pupil)
    featSig={vel,vel,acc,lick,rewardSignalC(sData,N),pos,pup};

    nCells=size(dff_all,1);
    for c=1:nCells
        of=famRes{c}; on=novRes{c};
        y=zscoreSafeC(double(dff_all(c,1:N))');
        if ~isempty(of)&&of.ok
            g=primaryCatC(of, signCorr(vel(famMask), y(famMask)));
            if g>=1
                [m,xc]=featureTuning(y, featSig{g}, famMask, nBins);
                GT.tunFam{g}=[GT.tunFam{g}; m]; GT.xcFam{g}=[GT.xcFam{g}; xc];
            end
        end
        if ~isempty(on)&&on.ok
            g=primaryCatC(on, signCorr(vel(novMask), y(novMask)));
            if g>=1
                [m,xc]=featureTuning(y, featSig{g}, novMask, nBins);
                GT.tunNov{g}=[GT.tunNov{g}; m]; GT.xcNov{g}=[GT.xcNov{g}; xc];
            end
        end
    end
end

% Mean dF/F of one cell binned by feature VALUE (quantile bins). Returns the
% per-bin dF/F (1 x nBins) and the per-bin feature-value centers (1 x nBins).
function [m, xc] = featureTuning(y, sig, mask, nBins)
    m=nan(1,nBins); xc=nan(1,nBins);
    if isempty(sig), return; end
    sig=sig(:); v=sig(mask); yy=y(mask);
    good=~isnan(v)&~isnan(yy); v=v(good); yy=yy(good);
    if numel(v)<nBins, return; end
    % quantile edges (robust to skwed feature distributions)
    edges=quantile(v, linspace(0,1,nBins+1));
    edges=unique(edges); if numel(edges)<3, return; end
    nb=numel(edges)-1;
    bi=discretize(v,edges);   % NaN (outside quantile range) -> excluded from all bins (bi==b is false)
    for b=1:nb
        sel=bi==b;
        if any(sel), m(b)=mean(yy(sel),'omitnan'); xc(b)=mean(v(sel),'omitnan'); end
    end
end

% =====================================================================
function plotCohortGroupedTuning(GT, cat_names)
    nG=GT.nG;
    figure('Name','Cohort feature tuning by winning feature','Color','w','Position',[60 60 820 1100]);
    sgtitle('Cohort feature tuning (pooled cells): each group''s \DeltaF/F vs ITS OWN feature', ...
        'FontWeight','bold','FontSize',12);
    c_fam=[0.2 0.5 0.8]; c_nov=[0.8 0.3 0.3];
    for gi=1:nG
        % Familiar
        subplot(nG,2,(gi-1)*2+1); hold on;
        drawFeatTuning(GT.xcFam{gi}, GT.tunFam{gi}, c_fam);
        ylabel({cat_names{gi},'\DeltaF/F (z)'},'FontWeight','bold','FontSize',8);
        if gi==1, title('Familiar'); end
        if gi==nG, xlabel(sprintf('feature value')); end
        text(0.02,0.86,sprintf('n=%d',size(GT.tunFam{gi},1)),'Units','normalized','FontSize',8);
        xlabel(GT.featName{gi},'FontSize',7);
        set(gca,'TickDir','out'); box off;
        % Novel
        subplot(nG,2,(gi-1)*2+2); hold on;
        drawFeatTuning(GT.xcNov{gi}, GT.tunNov{gi}, c_nov);
        if gi==1, title('Novel'); end
        text(0.02,0.86,sprintf('n=%d',size(GT.tunNov{gi},1)),'Units','normalized','FontSize',8);
        xlabel(GT.featName{gi},'FontSize',7);
        set(gca,'TickDir','out'); box off;
    end
end

% Plot pooled dF/F vs feature value: x = mean bin-center across cells, band =
% mean +/- SEM of dF/F across cells per bin.
function drawFeatTuning(XC, M, col)
    if isempty(M) || all(isnan(M(:))), text(0.4,0.5,'no cells','Units','normalized','Color',[.5 .5 .5]); return; end
    x = mean(XC,1,'omitnan');              % feature-value axis (mean center per bin)
    mu= mean(M,1,'omitnan'); n=sum(~isnan(M),1); sem=std(M,0,1,'omitnan')./sqrt(max(n,1));
    good=~isnan(x)&~isnan(mu);
    x=x(good); mu=mu(good); sem=sem(good);
    [x,ord]=sort(x); mu=mu(ord); sem=sem(ord);
    if numel(x)<2, plot(x,mu,'o','Color',col); return; end
    fill([x fliplr(x)],[mu+sem fliplr(mu-sem)],col,'FaceAlpha',0.25,'EdgeColor','none');
    plot(x,mu,'Color',col,'LineWidth',1.5);
end

function cat = primaryCatC(o, velSign)
    cat=0; sig=o.sig; dR2=o.dR2_unique; valid=find(sig & ~isnan(dR2));
    if isempty(valid), return; end
    [~,m]=max(dR2(valid)); pf=valid(m);
    switch pf
        case 1
            % velocity primary; if sign undefined (NaN), default to (+) so the
            % cell is still counted as velocity-primary rather than dropped.
            if isnan(velSign), cat=1; else, cat=(velSign>=0)*1+(velSign<0)*2; end
        case 2, cat=3; case 3, cat=4; case 4, cat=5; case 5, cat=6; case 6, cat=7;
    end
end

function progressTick(d, nCells, reset)
% Advances the progress dialog (called via DataQueue from parfor; runs on client).
    persistent done;
    if reset, done=0; return; end
    done=done+1;
    if ~isempty(d) && isvalid(d)
        d.Value=min(done/nCells,1);
        d.Message=sprintf('Fitting cells: %d/%d', min(done,nCells), nCells);
    end
end

function rew = rewardSignalC(sData, N)
    rew=zeros(N,1); vt=getSigC(sData,'rewardDs');
    if ~isempty(vt)&&any(double(vt)>0), vt=double(vt); rew(1:min(N,numel(vt)))=vt(1:min(N,numel(vt)))>0; return; end
    if isfield(sData,'daqdata')&&isfield(sData.daqdata,'frameIndex')&&isfield(sData.daqdata,'waterValve')
        vr=double(sData.daqdata.waterValve>0.5); ft=find(diff(sData.daqdata.frameIndex)==1)+1; lim=min(N,length(ft)-1);
        for i=1:lim, if any(vr(ft(i):ft(i+1)-1)), rew(i)=1; end, end
    end
end

function v = colvecC(x, N)
    if isempty(x), v=zeros(N,1); return; end
    x=double(x(:)); if numel(x)>=N, v=x(1:N); else, v=[x;zeros(N-numel(x),1)]; end
end

function s = getSigC(sData, name)
    if isfield(sData.behavior.signals,name), s=sData.behavior.signals.(name); else, s=[]; end
end

function PT = accumPositionTuning(PT, dff_all, sData, famRes, velTrace, pupilTrace, famMask, novMask, N)
    nBins=PT.nBins;
    pos=double(sData.behavior.signals.corridorPositionDs(1:N)); pos=pos(:);
    pmin=min(pos); pmax=max(pos); if pmax<=pmin, return; end
    posN=(pos-pmin)/(pmax-pmin); edges=linspace(0,1,nBins+1);
    binIdx=discretize(posN,edges);   % NaN (dropped/out-of-range) -> excluded from all bins (binIdx==b is false)
    trialIdx=double(sData.behavior.signals.trialIndexesDs(1:N)); trialIdx=trialIdx(:);

    vel=velTrace(:);
    timeA=double(sData.behavior.signals.timeDs(1:N)); FsA=1/mean(diff(timeA),'omitnan');
    acc=gradient(vel)*FsA;   % match GLM engine (cm/s^2)
    lick=colvecC(getSigC(sData,'lickIfreqDs'),N); rew=rewardSignalC(sData,N); pup=pupilTrace(:);
    featSig={vel,vel,acc,lick,rew,pos,pup};

    % per-session feature curves per group (mean across trials)
    for g=1:PT.nG
        PT.featFam{g}=[PT.featFam{g}; featBandMeanC(featSig{g},binIdx,trialIdx,famMask,nBins)];
        PT.featNov{g}=[PT.featNov{g}; featBandMeanC(featSig{g},binIdx,trialIdx,novMask,nBins)];
    end
    PT.rewFam{end+1}=rewardPositionsC(rew,posN,famMask);
    PT.rewNov{end+1}=rewardPositionsC(rew,posN,novMask);

    % per-cell dF/F position tuning, grouped by familiar winning feature
    nCells=size(dff_all,1);
    for c=1:nCells
        of=famRes{c}; if isempty(of)||~of.ok, continue; end
        y=zscoreSafeC(double(dff_all(c,1:N))');
        g=primaryCatC(of, signCorr(vel(famMask), y(famMask)));
        if g<1, continue; end
        PT.tunFam{g}=[PT.tunFam{g}; binMeanC(y,binIdx,famMask,nBins)];
        PT.tunNov{g}=[PT.tunNov{g}; binMeanC(y,binIdx,novMask,nBins)];
    end
end

function plotCohortPositionTuning(PT, cat_names)
    nG=PT.nG; nBins=PT.nBins; ctr=((1:nBins)-0.5)/nBins;
    figure('Color','w','Position',[60 40 820 1250]);
    sgtitle('Cohort position tuning (pooled cells) + mean feature; green = reward zones', ...
        'FontWeight','bold','FontSize',12);
    c_fam=[0.2 0.5 0.8]; c_nov=[0.8 0.3 0.3]; c_feat=[0.35 0.35 0.35];
    rewFam=meanRewLoc(PT.rewFam); rewNov=meanRewLoc(PT.rewNov);
    nRow=2*nG;
    for gi=1:nG
        rTop=(gi-1)*4;
        subplot(nRow,2,rTop+1); hold on; drawBandC(ctr,PT.tunFam{gi},c_fam);
        ylabel({cat_names{gi},'\DeltaF/F (z)'},'FontWeight','bold','FontSize',8);
        if gi==1, title('Familiar'); end
        text(0.02,0.85,sprintf('n=%d',size(PT.tunFam{gi},1)),'Units','normalized','FontSize',8);
        set(gca,'TickDir','out','XTickLabel',[]); box off; xlim([0 1]);
        subplot(nRow,2,rTop+2); hold on; drawBandC(ctr,PT.tunNov{gi},c_nov);
        if gi==1, title('Novel'); end
        text(0.02,0.85,sprintf('n=%d',size(PT.tunNov{gi},1)),'Units','normalized','FontSize',8);
        set(gca,'TickDir','out','XTickLabel',[]); box off; xlim([0 1]);
        subplot(nRow,2,rTop+3); hold on; drawBandC(ctr,PT.featFam{gi},c_feat); addRewardLinesC(rewFam);
        ylabel('feature','FontSize',8); if gi==nG, xlabel('Track position (norm.)'); end
        set(gca,'TickDir','out'); box off; xlim([0 1]);
        subplot(nRow,2,rTop+4); hold on; drawBandC(ctr,PT.featNov{gi},c_feat); addRewardLinesC(rewNov);
        if gi==nG, xlabel('Track position (norm.)'); end
        set(gca,'TickDir','out'); box off; xlim([0 1]);
    end
end

function m = binMeanC(y, binIdx, mask, nBins)
    m=nan(1,nBins);
    for b=1:nBins, sel=mask&(binIdx==b); if any(sel), m(b)=mean(y(sel),'omitnan'); end, end
end

function m = featBandMeanC(sig, binIdx, trialIdx, mask, nBins)
    m=nan(1,nBins); if isempty(sig), return; end
    trials=unique(trialIdx(mask & trialIdx>0)); if isempty(trials), return; end
    pt=nan(numel(trials),nBins);
    for tIx=1:numel(trials)
        tsel=mask&(trialIdx==trials(tIx));
        for b=1:nBins, s=tsel&(binIdx==b); if any(s), pt(tIx,b)=mean(sig(s),'omitnan'); end, end
    end
    m=mean(pt,1,'omitnan');
end

function drawBandC(x, M, col)
    if isempty(M)||all(isnan(M(:))), text(0.4,0.5,'no cells','Units','normalized','Color',[.5 .5 .5]); return; end
    mu=mean(M,1,'omitnan'); n=sum(~isnan(M),1); sem=std(M,0,1,'omitnan')./sqrt(max(n,1));
    good=~isnan(mu); xx=x(good); mu=mu(good); sem=sem(good);
    fill([xx fliplr(xx)],[mu+sem fliplr(mu-sem)],col,'FaceAlpha',0.25,'EdgeColor','none');
    plot(xx,mu,'Color',col,'LineWidth',1.5);
end

function addRewardLinesC(locs)
    if isempty(locs), return; end
    yl=ylim; for k=1:numel(locs), plot([locs(k) locs(k)],yl,'--','Color',[0.1 0.6 0.2],'LineWidth',1); end
    ylim(yl);
end

function locs = rewardPositionsC(rew, posN, mask)
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

function z = zscoreSafeC(y)
    s=std(y,'omitnan'); if s<1e-9, z=y-mean(y,'omitnan'); else, z=(y-mean(y,'omitnan'))/s; end
end

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

    ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2);
    ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
    trialIdx=sData.behavior.signals.trialIndexesDs(1:N); trialIdx=trialIdx(:);
    vtr=trialIdx>0 & trialIdx<=length(ctx);
    cctx=nan(N,1); cctx(vtr)=ctx(trialIdx(vtr)); cctx=fillmissing(cctx,'nearest');
    % Familiar id: prefer explicit field (robust if recording doesn't start familiar);
    % else first trial's context (correct for the standard fam->novel->fam protocol).
    if isfield(sData,'sessionInfo') && isfield(sData.sessionInfo,'famID')
        fam_id=sData.sessionInfo.famID;
    else
        fam_id=ctx(1);
    end
    uctx=unique(ctx); nov_id=uctx(uctx~=fam_id);
    famMask=(cctx==fam_id);   % all familiar frames (both blocks)
    if isempty(nov_id), novMask=false(N,1); else, novMask=(cctx==nov_id(1)); end
end

function closeIfValid(d), if isvalid(d), close(d); end, end

function s = ternStr(cond, a, b)
% tiny inline conditional for strings: ternStr(cond,a,b) -> a if cond else b
    if cond, s=a; else, s=b; end
end
