function vipGLM_paperFigure_dots()
% Copy of vipGLM_paperFigure.m (NatureComms Revision 2a) with one change:
% Fig. 2b shows the per-session values (n = 10 sessions) as dots on each bar,
% and the n= labels sit above the highest dot.
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_PAPERFIGURE  Build the VIP encoding-diversity figure (panels a-g).
%
% AIM: an encoding-model description of VIP response diversity that replaces the
% k-means cluster figure, answering the reviewer's request to quantify response
% diversity with a per-cell GLM (Turi et al. 2019 approach) rather than an
% arbitrary k=3 clustering.
%
% INPUT  (two file prompts):
%   STEP 1: the COHORT RESULTS file  -- ONE .mat (cohortGLM_*.mat, has res.cellTable)
%   STEP 2: the SESSION DATA files   -- ONE OR MORE sData .mat (multi-select)
% OUTPUT: one multi-panel MATLAB figure (panels a-g). Saves nothing.
% NEEDS:  vipGLM_drawExampleTrace.m on the path.
%
% LAYOUT (panels a-g):
%   a  per-cell unique dR2 heatmap   (viridis; features ranked by dR2,
%                                     cells block-sorted into a staircase)
%   b  mixed-selectivity histogram   (open bars, mean % across sessions
%                                     +/- SEM; inset: full-model CV R2)
%   c  example cell                  (observed vs held-out predicted dF/F + velocity)
%   d  unique dR2 per feature        (ranked high->low; velocity split Vel+ / Vel-;
%                                     merged-group rank-sum bracket)
%   e  2x2 tuning curves             (Vel+, Vel-, Pupil, Licking) + position-peak hist
%   f  mean+/-SEM dF/F vs trial      (small multiples: Pupil, Vel, Position, Licking)
%   g  upregulated cells in novel    (All + per feature; velocity split Vel+ / Vel-)
%
% CONVENTIONS (deliberate):
%   * Categories/feature stats from the FAMILIAR fit ONLY (full data).
%   * Novel appears only in f (trial traces) and g (upregulation), from raw dF/F.
%   * Cell order: features ranked by aggregate dR2 -> column order. Each cell
%     assigned to the block of its highest-dR2 feature; within block, sorted by
%     that feature's dR2 descending. This staircase order is shared by i and m.
%   * Colors: Okabe-Ito per feature (CVD-safe); viridis for dR2 magnitude (i);
%     diverging blue-white-red for signed dF/F (m). Arial, thick axes.
%
% INPUTS (prompted): STEP 1 cohort results .mat (res.cellTable);
%                    STEP 2 the session .mat files (multi-select).
% Serial (no parfor). NOTE: pupil panels require adequate pupil sampling -- the
% per-session detection fraction is printed; if low, pupil results are not
% interpretable. Friedman stats are read from the saved res.stats; if degenerate
% (NaN omnibus) the k panel prints a note instead of drawing brackets.

    FEATS  = {'Velocity','Acceleration','Licking','Reward','Position','Pupil'};
    ABBR   = {'Speed','Acc','Lick','Rew','Pos','Pupil'};
    % Okabe-Ito
    OI.vel=[0 114 178]/255;  OI.velL=[156 203 232]/255; OI.acc=[153 153 153]/255;
    OI.lick=[0 158 115]/255; OI.rew=[213 94 0]/255;     OI.pos=[230 159 0]/255;
    OI.pup=[204 121 167]/255;
    FCOL = [OI.vel; OI.acc; OI.lick; OI.rew; OI.pos; OI.pup];  % per-feature color

    % ---------------- load results ----------------
        uiwait(msgbox({'Select the COHORT RESULTS file.','','ONE file: the cohortGLM_*.mat saved by vipGLM_cohort (has res.cellTable).','','NOT a session/m8xxx file.'},'Cohort results (ONE file)','help'));
    [rf,rp]=uigetfile('*.mat','COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes);
    if isequal(rf,0), return; end
    R=load(fullfile(rp,rf));
    if ~isfield(R,'res')||~isfield(R.res,'cellTable'), error('No res.cellTable.'); end
    res=R.res; T=res.cellTable;
    Tf=T(T.context=="Familiar" & T.ok,:);
    hasVelSign = any(strcmpi(Tf.Properties.VariableNames,'velSign'));
    nCellF=height(Tf);

    % ---------------- per-cell dR2 matrix (familiar) ----------------
    DR2=nan(nCellF,6); SIG=false(nCellF,6);
    for f=1:6
        dc=['dR2_' FEATS{f}]; sc=['sig_' FEATS{f}];
        if any(strcmpi(Tf.Properties.VariableNames,dc)), DR2(:,f)=Tf.(dc); end
        if any(strcmpi(Tf.Properties.VariableNames,sc)), SIG(:,f)=Tf.(sc)>0; end
    end
    DR2d=DR2; DR2d(isnan(DR2d))=0; DR2d(DR2d<0)=0;   % display copy

    % ---------------- feature ranking (aggregate dR2) -> column order --------
    aggScore = sum(DR2d,1) + max(DR2d,[],1);          % total + peak
    [~,colOrder]=sort(aggScore,'descend');            % dominant feature first
    rankOfFeat=zeros(1,6); rankOfFeat(colOrder)=1:6;

    % ---------------- staircase cell order (shared by i and m) ---------------
    % primary feature = argmax dR2 AMONG SIGNIFICANT features (matches cohort
    % primIdx definition). Cells with no significant feature fall back to the
    % global argmax purely so they still get a heatmap row; they are excluded
    % from l/m/n by the isSigPrimary guard anyway.
    [primV,primF]=max(DR2d,[],2);                 % global argmax (fallback / heatmap sort)
    for gi=1:nCellF
        sigf=find(SIG(gi,:));
        if ~isempty(sigf)
            [pv,mm]=max(DR2d(gi,sigf));
            primF(gi)=sigf(mm); primV(gi)=pv;     % argmax within significant set
        end
    end
    sortKey = rankOfFeat(primF)' * 1e6 - primV;        % block first, then dR2 desc
    [~,cellOrder]=sort(sortKey,'ascend');

    % ---------------- velocity dR2 split by sign: accumulated LIVE in loop ----
    % (collect each Vel-primary cell's dR2_Velocity into +/- pools using the
    %  live signCorr sign, so k matches l/m/n grouping; mean+SEM computed after)
    velPpool=[]; velMpool=[];

    % ---------------- per-feature mean dR2 (sig cells) for k ------------------
    fMean=nan(1,6); fSEM=nan(1,6);
    for f=1:6
        v=DR2(SIG(:,f),f); v=v(~isnan(v));
        if ~isempty(v), fMean(f)=mean(v); fSEM(f)=std(v)/sqrt(numel(v)); end
    end

    % ---------------- mixed selectivity: per-session % (for mean/SEM) --------
    sess=unique(string(Tf.session),'stable'); nSess=numel(sess);
    nSig=sum(SIG,2);                              % per cell
    muxBins=0:4;                                  % collapse >=4 into 4
    nSigC=min(nSig,4);
    muxPerSess=nan(nSess,numel(muxBins));
    for s=1:nSess
        sel=string(Tf.session)==sess(s); ns=sum(sel);
        if ns==0, continue; end
        for b=1:numel(muxBins)
            muxPerSess(s,b)=100*sum(nSigC(sel)==muxBins(b))/ns;
        end
    end
    muxMean=mean(muxPerSess,1,'omitnan');
    muxSEM =std(muxPerSess,0,1,'omitnan')./sqrt(sum(~isnan(muxPerSess),1));

    % ---------------- Friedman / Holm from saved stats (familiar) ------------
    frP=NaN; pairIdx=[]; pairHolm=[]; statNote='';
    if isfield(res,'stats') && isfield(res.stats,'dR2_fam')
        sf=res.stats.dR2_fam;
        if isfield(sf,'friedman_p'), frP=sf.friedman_p; end
        if isfield(sf,'pair_idx'),  pairIdx=sf.pair_idx; end
        if isfield(sf,'pair_p_holm'),pairHolm=sf.pair_p_holm; end
        if isfield(sf,'note'),       statNote=char(string(sf.note)); end
    end

    % ===================== SESSIONS: raster, tuning, upregulation ============
        uiwait(msgbox({'Select the SESSION DATA files.','','ONE OR MORE m8xxx session .mat files (hold Cmd/Shift for several).'},'Session data (one or more)','help'));
    [sfn,sp]=uigetfile('*.mat','SESSION DATA files (ONE OR MORE m8xxx .mat; multi-select)',dfSess,'MultiSelect','on');
    if isequal(sfn,0), return; end
    if ischar(sfn), sfn={sfn}; end

    nGrid=110; nbin=20;
    rastTrace=nan(nCellF,nGrid);                   % per cell (cellTable order)
    tun=struct('velP',[],'velM',[],'pup',[],'lick',[]);
    tunX=struct('velP',[],'velM',[],'pup',[],'lick',[]);   % per-bin feature-value centers
    mTrace=struct('velP',[],'velM',[],'lick',[],'pos',[],'pup',[]);
    posPeak=[];   % peak track position (cm) for each Position-primary cell
    upInc=zeros(1,6); upTot=zeros(1,6);            % per feature (Vel combined here)
    upIncV=[0 0]; upTotV=[0 0];                    % [Vel+ ; Vel-]
    upTotAll=0; upIncAll=0;                         % all cells, unguarded (for the 'All' bar)
    pupDet=[];

    keySess=string(Tf.session); keyID=Tf.cellID;
    % ----- example cell shown in the panel-k strip (observed vs held-out GLM) -----
    exCfg=struct('sessTag','m8130-20231024-01-K5', ...   % substring matching the session filename
                 'cellID',1, 'tWin',[7 37], 'lpfSec',0.25);
    exSData=[];
    for si=1:numel(sfn)
        fn=sfn{si}; S=load(fullfile(sp,fn));
        if ~isfield(S,'sData'), warning('skip %s',fn); continue; end
        sData=S.sData;
        if contains(fn, exCfg.sessTag), exSData=sData; end   % stash for panel-k example trace
        if ~(isfield(sData.imdata,'roiSignals')&&numel(sData.imdata.roiSignals)>=2)
            warning('skip %s (no roiSignals(2))',fn); continue; end
        dff=sData.imdata.roiSignals(2).dff; nCells=size(dff,1); nF=size(dff,2);
        ti=double(sData.behavior.signals.trialIndexesDs); ti=ti(:); N=min(nF,numel(ti)); ti=ti(1:N);
        vel=double(sData.behavior.signals.velocityDs(1:N)); vel=vel(:);
        if numel(vel)>2, vel(1:2)=vel(3); end   % match cohort quickBehavior: drop startup spike
        pos=double(sData.behavior.signals.corridorPositionDs(1:N)); pos=pos(:);
        lick=colvec(getSig(sData,'lickIfreqDs'),N);
        pup =pupilTrace(sData,N);
        ctx=round(fillmissing(mode(sData.behavior.trialMatrices.contextMatrix,2),'nearest')); ctx=ctx(:);
        fam_id=ctx(1); uc=unique(ctx); nov_id=uc(uc~=fam_id);
        vtr=ti>0&ti<=length(ctx); cc=nan(N,1); cc(vtr)=ctx(ti(vtr)); cc=fillmissing(cc,'nearest');
        famMask=(cc==fam_id);
        if isempty(nov_id), warning('skip %s (no novel)',fn); continue; end
        pupDet(end+1)=mean(~isnan(pup(famMask))); %#ok<AGROW>

        allTr=unique(ti(ti>0));
        isFam=arrayfun(@(t) mode(cc(ti==t))==fam_id, allTr);
        novTr=allTr(~isFam); if isempty(novTr), continue; end
        firstNov=novTr(1);
        initFam=allTr(isFam & allTr<firstNov);
        lastFam=initFam(max(1,end-9):end); firstNov10=novTr(1:min(10,end));

        rows=find(keySess==string(fn));
        for rr=1:numel(rows)
            gi=rows(rr); cid=keyID(gi); if cid<1||cid>nCells, continue; end
            y=double(dff(cid,1:N))'; yz=zscoreSafe(y);
            % velocity sign from stored cellTable velSign (default + if absent/NaN)
            vsgn=1; if hasVelSign && ~isnan(Tf.velSign(gi)), vsgn=sign(Tf.velSign(gi)); end
            if vsgn==0, vsgn=1; end
            % raster trace (per-trial mean, resampled), placed at cellTable row gi
            tv=nan(1,numel(allTr));
            for tk=1:numel(allTr), s2=ti==allTr(tk); if any(s2), tv(tk)=mean(yz(s2),'omitnan'); end, end
            rastTrace(gi,:)=resampleToGrid(double(allTr(:)'),tv,firstNov,nGrid);
            pf=primF(gi);
            isSigPrimary = SIG(gi,pf);   % dilution guard: cell genuinely significant for its top feature

            % novelty upregulation (Mate's criterion; adjacent-window drift control)
            fM=arrayfun(@(t) mean(y(ti==t),'omitnan'), lastFam); fM=fM(~isnan(fM));
            nM=arrayfun(@(t) mean(y(ti==t),'omitnan'), firstNov10); nM=nM(~isnan(nM));
            if numel(fM)>=3 && numel(nM)>=3
                np=min(numel(fM),numel(nM)); inc=false;
                fM_adj=fM(end-np+1:end);   % anchor familiar to END (closest to novel onset)
                nM_adj=nM(1:np);
                try [~,pp]=ttest(nM_adj,fM_adj,'Tail','right'); inc=pp<0.05; catch, inc=mean(nM_adj)>mean(fM_adj); end
                upTotAll=upTotAll+1; if inc, upIncAll=upIncAll+1; end   % ALL cells (unguarded)
                if isSigPrimary
                    upTot(pf)=upTot(pf)+1; if inc, upInc(pf)=upInc(pf)+1; end
                    if pf==1
                        k=(vsgn<0)+1;  % 1=Vel+ , 2=Vel-
                        upTotV(k)=upTotV(k)+1; if inc, upIncV(k)=upIncV(k)+1; end
                    end
                end
            end

            % everything below summarises a cell by its PRIMARY feature -> require significance
            if ~isSigPrimary, continue; end

            % peak track position for cells whose PRIMARY feature is position (familiar)
            if pf==5
                [pm,px]=featTuneBins(yz, pos, famMask, 0, 250, 40);  % finer spatial bins
                if any(~isnan(pm)), [~,ix]=max(pm); posPeak(end+1)=px(ix); end %#ok<AGROW>
            end

            % familiar tuning by native variable (QUANTILE bins + centers, matching cohort)
            if pf==1
                if vsgn<0, velMpool(end+1)=DR2d(gi,1); else, velPpool(end+1)=DR2d(gi,1); end %#ok<AGROW>
                [tm,tx]=featTuneQuantile(yz, vel, famMask, nbin);
                if vsgn<0, tun.velM=[tun.velM; tm]; tunX.velM=[tunX.velM; tx];
                else,      tun.velP=[tun.velP; tm]; tunX.velP=[tunX.velP; tx]; end
            elseif pf==3
                [tm,tx]=featTuneQuantile(yz, lick, famMask, nbin);
                tun.lick=[tun.lick; tm]; tunX.lick=[tunX.lick; tx];
            elseif pf==6
                [tm,tx]=featTuneQuantile(yz, pup, famMask, nbin);
                tun.pup=[tun.pup; tm]; tunX.pup=[tunX.pup; tx];
            end
            % per-category trial-trace for panel m
            if pf==1
                if vsgn<0, mTrace.velM=[mTrace.velM; rastTrace(gi,:)];
                else,      mTrace.velP=[mTrace.velP; rastTrace(gi,:)]; end
            elseif pf==3, mTrace.lick=[mTrace.lick; rastTrace(gi,:)];
            elseif pf==5, mTrace.pos =[mTrace.pos;  rastTrace(gi,:)];
            elseif pf==6, mTrace.pup =[mTrace.pup;  rastTrace(gi,:)];
            end
        end
    end

    % velocity dR2 means + combined SEM (over all Vel-primary cells, both signs)
    velP_dR2=0; velM_dR2=0;
    if ~isempty(velPpool), velP_dR2=mean(velPpool,'omitnan'); end
    if ~isempty(velMpool), velM_dR2=mean(velMpool,'omitnan'); end
    velAll=[velPpool velMpool];
    velSEM=0; if numel(velAll)>1, velSEM=std(velAll,'omitnan')/sqrt(numel(velAll)); end

    % merged-group comparison: (Pupil + Velocity) vs (Lick + Position + Acceleration + Reward), no correction
    % pool each cell's dR2 in the relevant feature among significant cells
    poolPV=[DR2(SIG(:,6),6); DR2(SIG(:,1),1)];                              % Pupil + Velocity
    poolPL=[DR2(SIG(:,3),3); DR2(SIG(:,5),5); DR2(SIG(:,2),2); DR2(SIG(:,4),4)]; % Lick+Pos+Acc+Rew
    poolPV=poolPV(~isnan(poolPV)); poolPL=poolPL(~isnan(poolPL));
    pPVvsPL=NaN;
    if numel(poolPV)>1 && numel(poolPL)>1
        try, pPVvsPL=ranksum(poolPV,poolPL); catch, pPVvsPL=NaN; end
    end

    % ===================== FIGURE =====================
    set(0,'DefaultAxesFontName','Arial','DefaultTextFontName','Arial');
    fig=figure('Color','w','Position',[30 30 1300 760],'Renderer','painters');

    % explicit axes positions so i/j and m/n share widths exactly
    % columns: left x=0.06 w=0.20 ; mid x=0.36 w=0.26 ; right x=0.72 w=0.24
    Lx=0.060; Lw=0.200; Mx=0.345; Mw=0.275; Rx=0.715; Rw=0.235;
    topY=0.560; topH=0.380; botY=0.085; botH=0.360;
    LW=1.4; FS=10; FSlab=11;

    % ---- i: heatmap (viridis), ranked columns, staircase rows ----
    axI=axes('Position',[Lx topY Lw topH]);
    Mi=DR2d(cellOrder, colOrder);
    imagesc(axI, Mi, [0 0.25]); colormap(axI, viridisMap());
    set(axI,'YDir','reverse','XTick',1:6,'XTickLabel',ABBR(colOrder),'YTick',[1 nCellF], ...
        'TickDir','out','LineWidth',LW,'FontSize',FS,'XTickLabelRotation',30);
    ylabel(axI,sprintf('Cells (n=%d)',nCellF),'FontSize',FSlab);
    title(axI,'a','FontSize',13,'FontWeight','bold','HorizontalAlignment','left','Units','normalized','Position',[-0.18 1.02 0]);
    cI=colorbar(axI,'Location','southoutside'); cI.Label.String='unique \DeltaR^2'; cI.LineWidth=LW; cI.FontSize=FS;

    % ---- j: mixed selectivity, black open bars, mean +/- SEM (gap after 0) ----
    % j occupies the full lower-left slot. Bars on the left ~60%, CV inset on the right.
    jW=Lw*0.60;
    jStretch=0.050;                               % taller main bars, grows upward
    axJ=axes('Position',[Lx botY jW botH+jStretch]); hold(axJ,'on');
    jx=muxBins + [0 0.4 0.4 0.4 0.4];   % shift 1..4 right to separate the 0-bar
    for b=1:numel(jx)
        if b==1
            bar(axJ, jx(b), muxMean(b), 0.6, 'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',LW);
        else
            bar(axJ, jx(b), muxMean(b), 0.6, 'FaceColor','none','EdgeColor','k','LineWidth',LW);
        end
    end
    rs=rng; rng(1);                              % fixed horizontal jitter of the session dots
    for b=1:numel(jx)
        v=muxPerSess(:,b); v=v(~isnan(v));
        scatter(axJ, jx(b)+(rand(numel(v),1)-0.5)*0.32, v, 9, [0.35 0.35 0.35], 'filled');
    end
    rng(rs);
    errorbar(axJ, jx, muxMean, muxSEM, 'k','LineStyle','none','LineWidth',LW,'CapSize',4);
    set(axJ,'TickDir','out','LineWidth',LW,'FontSize',FS,'XTick',jx,'XTickLabel',{'0','1','2','3','4'});
    xlabel(axJ,'Nr. of significant features per VIP cell','FontSize',FSlab); ylabel(axJ,'% of cells (\pm SEM)','FontSize',FSlab);
    box(axJ,'off'); xlim(axJ,[-0.6 5.0]);
    jN=arrayfun(@(b) sum(nSigC==b), muxBins);
    for b=1:numel(jx)
        text(axJ,jx(b),max([muxMean(b)+muxSEM(b); muxPerSess(:,b)])+3,sprintf('n=%d',jN(b)),'HorizontalAlignment','center','FontSize',FS);
    end
    title(axJ,'b','FontSize',13,'FontWeight','bold','HorizontalAlignment','left','Units','normalized','Position',[-0.30 1.02 0]);
    % CV R^2 inset: to the RIGHT of j's bars, same slot
    insW=Lw*0.30; insH=botH*0.46;
    insX=Lx+Lw-insW; insY=botY+botH-insH-botH*0.02;
    axJins=axes('Position',[insX insY insW insH]); hold(axJins,'on');
    r2=[];
    if any(strcmpi(Tf.Properties.VariableNames,'R2_full')), r2=Tf.R2_full; r2=r2(~isnan(r2)); end
    if ~isempty(r2)
        edges=-0.1:0.05:0.6; histogram(axJins,r2,edges,'FaceColor',[0.5 0.5 0.5],'EdgeColor','k','LineWidth',0.3);
    end
    set(axJins,'TickDir','out','LineWidth',1,'FontSize',FS-1); box(axJins,'off');
    xlabel(axJins,'Full-model CV R^2','FontSize',FS); ylabel(axJins,'# cells','FontSize',FS);
    xlim(axJins,[-0.1 0.6]); xticks(axJins,0:0.2:0.6);

    % ---- k: per-feature dR2 ranked, UNFILLED bars, velocity split via divider ----
    % k is compressed into the lower part of the middle-top slot; the top ~42% of
    % the slot (full width) is left EMPTY for a future actual-vs-predicted dF/F trace.
    kEmptyFrac=0.57;                              % top fraction reserved (taller k trace strip)
    kBandH=topH*(1-kEmptyFrac-0.06);              % height of the k bar band
    axK=axes('Position',[Mx topY Mw kBandH]); hold(axK,'on');
    kfeat=colOrder;                                  % ranked order
    yu=zeros(6,1); kN=zeros(6,1); xPos=nan(1,6);
    for c=1:6
        f=kfeat(c);
        if f==1
            bar(axK,c,velP_dR2+velM_dR2,0.6,'FaceColor','none','EdgeColor','k','LineWidth',LW);
            plot(axK,[c-0.3 c+0.3],[velP_dR2 velP_dR2],'k-','LineWidth',0.75);
            errorbar(axK,c,velP_dR2+velM_dR2,velSEM,'k','LineStyle','none','LineWidth',LW,'CapSize',4);
            % pos/neg labels inside the two velocity segments
            text(axK,c,velP_dR2/2,'pos','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FS-2);
            text(axK,c,velP_dR2+velM_dR2/2,'neg','Rotation',90,'HorizontalAlignment','center','VerticalAlignment','middle','FontSize',FS-2);
            yu(c)=velP_dR2+velM_dR2+velSEM; kN(c)=sum(SIG(:,1));
        else
            bar(axK,c,fMean(f),0.6,'FaceColor','none','EdgeColor','k','LineWidth',LW);
            if isfinite(fSEM(f)) && fSEM(f)>0
                errorbar(axK,c,fMean(f),fSEM(f),'k','LineStyle','none','LineWidth',LW,'CapSize',4);
            end
            yu(c)=fMean(f)+max(fSEM(f),0);
            if ~isfinite(yu(c)), yu(c)=max(fMean(f),0); end
            kN(c)=sum(SIG(:,f));
        end
        xPos(f)=c;   % map feature index -> x location
    end
    set(axK,'XTick',1:6,'XTickLabel',ABBR(kfeat),'TickDir','out','LineWidth',LW,'FontSize',FS,'XTickLabelRotation',30);
    ylabel(axK,'Unique \DeltaR^2','FontSize',FSlab); box(axK,'off'); xlim(axK,[0.4 6.6]); ylim(axK,[0 max(yu)*1.62]);
    for c=1:6, text(axK,c,yu(c)+max(yu)*0.10,sprintf('n=%d',kN(c)),'HorizontalAlignment','center','FontSize',FS); end
    % merged-group comparison bracket: (Pup+Vel) vs (Lick+Pos+Acc+Rew)
    if isfinite(pPVvsPL)
        ytop=max(yu);
        yb=ytop*1.22;           % horizontal bar height
        prong=ytop*0.05;        % short vertical prongs down onto each bar
        leftX =[xPos(6) xPos(1)];                 % Pupil, Velocity
        rightX=[xPos(3) xPos(5) xPos(2) xPos(4)]; % Lick, Pos, Acc, Rew
        % left group: prongs + spanning bar
        for xx=leftX,  plot(axK,[xx xx],[yb-prong yb],'k-','LineWidth',1); end
        plot(axK,[min(leftX) max(leftX)],[yb yb],'k-','LineWidth',1);  xL=mean(leftX);
        % right group: prongs + spanning bar
        for xx=rightX, plot(axK,[xx xx],[yb-prong yb],'k-','LineWidth',1); end
        plot(axK,[min(rightX) max(rightX)],[yb yb],'k-','LineWidth',1); xR=mean(rightX);
        % connector between the two group midpoints, one notch higher
        yc=yb+ytop*0.10;
        plot(axK,[xL xL],[yb yc],'k-','LineWidth',1);
        plot(axK,[xR xR],[yb yc],'k-','LineWidth',1);
        plot(axK,[xL xR],[yc yc],'k-','LineWidth',1);
        text(axK,mean([xL xR]),yc+ytop*0.12,sprintf('p=%.3g',pPVvsPL),'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',FS-1);
    end
    title(axK,'d','FontSize',13,'FontWeight','bold','HorizontalAlignment','left','Units','normalized','Position',[-0.14 1.02 0]);

    % ---- embedded example-cell trace in the reserved strip ABOVE k ----
    % strip spans the top (1-kEmptyFrac-...) of the middle-top slot; split into
    % a dF/F panel (upper ~65%) and a velocity panel (lower ~35%).
    stripY0=topY+kBandH+0.075;                         % bottom of the strip (gap above l)
    stripH =topY+topH-stripY0;                         % available height up to slot top
    if stripH>0.04 && ~isempty(exSData)
        velFrac=0.32; vGap=0.014;
        axVel=axes('Position',[Mx stripY0 Mw stripH*velFrac]);
        axTrc=axes('Position',[Mx stripY0+stripH*velFrac+vGap Mw stripH*(1-velFrac)-vGap]);
        okEx=vipGLM_drawExampleTrace(axTrc,axVel,exSData,exCfg.cellID,exCfg.tWin,exCfg.lpfSec,FSlab);
        if okEx
            % panel letter as a figure-level annotation, pinned to the SAME absolute
            % x as l and m (Mx + (-0.14)*Mw) so the j/k/l/m/n/o letters form a column
            kLetX = Mx - 0.14*Mw;
            kLetY = stripY0+stripH-0.004;     % just above the trace panel top
            annotation(fig,'textbox',[kLetX kLetY 0.03 0.03],'String','c', ...
                'FontSize',13,'FontWeight','bold','EdgeColor','none', ...
                'HorizontalAlignment','left','VerticalAlignment','middle');
        end
    end

    % ---- l: 2x2 tuning (top) + wide Position peak histogram (bottom) ----
    % top 2x2 occupies upper ~62% of the l slot; position histogram spans bottom
    lTopH=botH*0.60; lBotH=botH*0.26; lGap=0.05;
    lTopY=botY+lBotH+lGap;
    tw=(Mw-0.03)/2; rowGap=0.075; th=(lTopH-rowGap)/2;
    drawTuning(axes('Position',[Mx         lTopY+th+rowGap tw th]), tun.velP, mean(tunX.velP,1,'omitnan'), 'Speed (cm/s)',sprintf('Speed(+) (n=%d)',size(tun.velP,1)),LW,FS,OI.vel);
    drawTuning(axes('Position',[Mx+tw+0.03 lTopY+th+rowGap tw th]), tun.velM, mean(tunX.velM,1,'omitnan'), 'Speed (cm/s)',sprintf('Speed(-) (n=%d)',size(tun.velM,1)),LW,FS,OI.vel);
    drawTuning(axes('Position',[Mx         lTopY         tw th]), tun.pup,  mean(tunX.pup,1,'omitnan'),  'Pupil diam. (a.u.)',sprintf('Pupil (n=%d)',size(tun.pup,1)),LW,FS,OI.pup);
    axLick=axes('Position',[Mx+tw+0.03 lTopY         tw th]);
    drawTuning(axLick, tun.lick, mean(tunX.lick,1,'omitnan'),  'Lick rate (Hz)',sprintf('Licking (n=%d)',size(tun.lick,1)),LW,FS,OI.lick);
    xlim(axLick,[0 8]);   % cap at mouse max lick rate
    % wide Position histogram across the bottom
    axPos=axes('Position',[Mx botY Mw lBotH]); hold(axPos,'on');
    if ~isempty(posPeak)
        edges=0:15:255; histogram(axPos,posPeak,edges,'FaceColor',OI.pos,'EdgeColor','k','LineWidth',0.3);
    end
    rewardPos=[120 210];   % reward zone locations (cm) on the 250 cm track; edit if needed
    yl=ylim(axPos);
    for rp=rewardPos, plot(axPos,[rp rp],yl,'k--','LineWidth',1); end
    ylim(axPos,yl);
    set(axPos,'TickDir','out','LineWidth',LW,'FontSize',FS-2); box(axPos,'off'); xlim(axPos,[0 250]);
    xlabel(axPos,'Peak track position (cm)','FontSize',FS-1); ylabel(axPos,'# position cells','FontSize',FS-1);
    title(axPos,sprintf('Position (n=%d)',numel(posPeak)),'FontSize',FS,'FontWeight','bold','Color','k','HorizontalAlignment','left','Units','normalized','Position',[0 0.88 0]);
    axL=axes('Position',[Mx botY Mw botH],'Color','none','XColor','none','YColor','none');
    title(axL,'e','FontSize',13,'FontWeight','bold','HorizontalAlignment','left','Units','normalized','Position',[-0.14 1.04 0]);

    % ---- m: small multiples, mean+/-SEM dF/F vs trial (F->N->F), one per category ----
    % velocity combined (Vel+ and Vel- pooled); lick now occupies the freed panel
    Mlist={ mTrace.pup, [mTrace.velP; mTrace.velM], mTrace.pos, mTrace.lick };
    mcolS=[OI.pup; OI.vel; OI.pos; OI.lick];
    mlabS={'Pupil','Speed','Position','Licking'};
    xt=1:nGrid;
    % shared y-range across the four (computed from smoothed means+sem)
    ylo=inf; yhi=-inf; MU=cell(1,4); SE=cell(1,4); NN=zeros(1,4);
    for q=1:4
        Mq=Mlist{q};
        if isempty(Mq)||all(isnan(Mq(:))), MU{q}=[]; continue; end
        mu=smoothNan(mean(Mq,1,'omitnan'),7); ns=sum(~isnan(Mq),1);
        sem=smoothNan(std(Mq,0,1,'omitnan')./sqrt(max(ns,1)),7);
        MU{q}=mu; SE{q}=sem; NN(q)=size(Mq,1);
        ylo=min(ylo,min(mu-sem)); yhi=max(yhi,max(mu+sem));
    end
    if ~isfinite(ylo), ylo=-0.5; yhi=0.5; end
    pad=0.05*(yhi-ylo); ylo=ylo-pad; yhi=yhi+pad;
    mtw=(Rw-0.03)/2; mth=(topH-0.05)/2;
    mpos=[Rx           topY+mth+0.05 mtw mth;
          Rx+mtw+0.03  topY+mth+0.05 mtw mth;
          Rx           topY          mtw mth;
          Rx+mtw+0.03  topY          mtw mth];
    for q=1:4
        axq=axes('Position',mpos(q,:)); hold(axq,'on'); %#ok<LAXES>
        if ~isempty(MU{q})
            good=~isnan(MU{q}); xx=xt(good); mu=MU{q}(good); sem=SE{q}(good);
            fill(axq,[xx fliplr(xx)],[mu+sem fliplr(mu-sem)],mcolS(q,:),'FaceAlpha',0.20,'EdgeColor','none');
            plot(axq,xx,mu,'Color',mcolS(q,:),'LineWidth',1.6);
            for xv=[30 90]/110*nGrid, plot(axq,[xv xv],[ylo yhi],'k--','LineWidth',0.75); end
        end
        ylim(axq,[ylo yhi]); xlim(axq,[1 nGrid]);
        set(axq,'TickDir','out','LineWidth',LW,'FontSize',FS-2,'XTick',[1 30 90 110]/110*nGrid,'XTickLabel',{'0','30','90','110'});
        title(axq,mlabS{q},'FontSize',FS,'FontWeight','bold','Color','k','HorizontalAlignment','left','Units','normalized','Position',[0 1.0 0]);
        if q==1||q==3, ylabel(axq,'\DeltaF/F (z)','FontSize',FS-1); end
        if q>=3, xlabel(axq,'Trial','FontSize',FS-1); end
    end
    axMlbl=axes('Position',[Rx topY Rw topH],'Color','none','XColor','none','YColor','none');
    title(axMlbl,'f','FontSize',13,'FontWeight','bold','HorizontalAlignment','left','Units','normalized','Position',[-0.15 1.04 0]);

    % ---- n: % upregulated cells in novel; first bar = all VIP cells, then per feature ----
    %        UNFILLED bars; SEM = binomial sqrt(p(1-p)/n); n shown above bars
    axN=axes('Position',[Rx botY Rw botH]); hold(axN,'on');
    % bar 1: all VIP cells (pooled across every cell with a valid novelty window)
    nAll=max(upTotAll,1); pAll=upIncAll/nAll;
    pctAll=100*pAll; seAll=100*sqrt(pAll*(1-pAll)/nAll);
    bar(axN,1,pctAll,0.6,'FaceColor',[0.8 0.8 0.8],'EdgeColor','k','LineWidth',LW);
    errorbar(axN,1,pctAll,seAll,'k','LineStyle','none','LineWidth',LW,'CapSize',4);
    text(axN,1,min(pctAll+seAll+4,99),sprintf('%d/%d',upIncAll,nAll),'HorizontalAlignment','center','FontSize',FS);
    % bars 2..: per feature in ranked order, shifted right to separate from 'All'.
    % Velocity is split into two bars (Vel+ and Vel-) using the signed counts.
    xpos=[]; labs={}; x=2+0.4;               % running x position for feature bars
    for c=1:6
        f=kfeat(c);
        if f==1
            % two velocity bars: Vel+ then Vel-
            velLab={'Speed(+)','Speed(-)'};
            for s=1:2
                ntot=max(upTotV(s),1); nup=upIncV(s);
                p=nup/ntot; pct=100*p; se=100*sqrt(p*(1-p)/ntot);
                bar(axN,x,pct,0.6,'FaceColor','none','EdgeColor','k','LineWidth',LW);
                errorbar(axN,x,pct,se,'k','LineStyle','none','LineWidth',LW,'CapSize',4);
                text(axN,x,min(pct+se+4,99),sprintf('%d/%d',nup,ntot),'HorizontalAlignment','center','FontSize',FS);
                xpos(end+1)=x; labs{end+1}=velLab{s}; %#ok<AGROW>
                x=x+1;
            end
        else
            ntot=max(upTot(f),1); nup=upInc(f);
            p=nup/ntot; pct=100*p; se=100*sqrt(p*(1-p)/ntot);
            bar(axN,x,pct,0.6,'FaceColor','none','EdgeColor','k','LineWidth',LW);
            errorbar(axN,x,pct,se,'k','LineStyle','none','LineWidth',LW,'CapSize',4);
            text(axN,x,min(pct+se+4,99),sprintf('%d/%d',nup,ntot),'HorizontalAlignment','center','FontSize',FS);
            xpos(end+1)=x; labs{end+1}=ABBR{f}; %#ok<AGROW>
            x=x+1;
        end
    end
    set(axN,'XTick',[1 xpos],'XTickLabel',[{'All'} labs],'TickDir','out','LineWidth',LW,'FontSize',FS,'XTickLabelRotation',30);
    ylabel(axN,'% upregulated VIP cells in novel','FontSize',FSlab); box(axN,'off'); xlim(axN,[0.4 max(xpos)+1]); ylim(axN,[0 100]);
    title(axN,'g','FontSize',13,'FontWeight','bold','HorizontalAlignment','left','Units','normalized','Position',[-0.15 1.02 0]);

    set(fig,'Name','VIP encoding diversity (a-g)','NumberTitle','off');
end

% ===================== helpers =====================
function drawTuning(ax, M, xc, xlab, catlab, LW, FS, col)
    if nargin<8 || isempty(col), col=[0 0 0]; end
    hold(ax,'on');
    if ~isempty(M) && ~all(isnan(M(:)))
        mu=mean(M,1,'omitnan'); n=sum(~isnan(M),1); sem=std(M,0,1,'omitnan')./sqrt(max(n,1));
        good=~isnan(mu)&~isnan(xc); xx=xc(good); mu=mu(good); sem=sem(good);
        [xx,o]=sort(xx); mu=mu(o); sem=sem(o);
        if numel(xx)>=2
            fill(ax,[xx fliplr(xx)],[mu+sem fliplr(mu-sem)],col,'FaceAlpha',0.20,'EdgeColor','none');
            plot(ax,xx,mu,'Color',col,'LineWidth',2);
        end
    else
        text(ax,0.5,0.5,'no cells','Units','normalized','HorizontalAlignment','center','Color',[0.6 0.6 0.6],'FontSize',8);
    end
    set(ax,'TickDir','out','LineWidth',LW,'FontSize',FS-2); box(ax,'off');
    xlabel(ax,xlab,'FontSize',FS-1); ylabel(ax,'\DeltaF/F (z)','FontSize',FS-1);
    title(ax,catlab,'FontSize',FS,'FontWeight','bold','Color','k','HorizontalAlignment','left','Units','normalized','Position',[0 1.0 0]);
end

function [m, xc] = featTuneBins(y, v, mask, lo, hi, nb)
    m=nan(1,nb); xc=nan(1,nb);
    v=v(mask); yy=y(mask); good=~isnan(v)&~isnan(yy); v=v(good); yy=yy(good);
    if isempty(v), return; end
    edges=linspace(lo,hi,nb+1); bi=discretize(v,edges);
    for b=1:nb
        sel=bi==b;
        if any(sel), m(b)=mean(yy(sel),'omitnan'); end
        xc(b)=(edges(b)+edges(b+1))/2;
    end
end

function [m, xc] = featTuneQuantile(y, sig, mask, nBins)
    m=nan(1,nBins); xc=nan(1,nBins);
    if isempty(sig), return; end
    sig=sig(:); v=sig(mask); yy=y(mask);
    good=~isnan(v)&~isnan(yy); v=v(good); yy=yy(good);
    if numel(v)<nBins, return; end
    edges=quantile(v, linspace(0,1,nBins+1)); edges=unique(edges);
    if numel(edges)<3, return; end
    nb=numel(edges)-1; bi=discretize(v,edges); bi(isnan(bi))=nb;
    for b=1:nb
        sel=bi==b;
        if any(sel), m(b)=mean(yy(sel),'omitnan'); xc(b)=mean(v(sel),'omitnan'); end
    end
end

function ys = smoothNan(y, w)
    ys=y; n=numel(y); h=floor(w/2);
    for i=1:n
        lo=max(1,i-h); hi=min(n,i+h); seg=y(lo:hi);
        if any(~isnan(seg)), ys(i)=mean(seg,'omitnan'); end
    end
end

function grid = resampleToGrid(tn, tv, novFirstTrial, nGrid)
    grid=nan(1,nGrid); good=~isnan(tv); if nnz(good)<3, return; end
    tn=tn(good); tv=tv(good); t0=min(tn); tNov=novFirstTrial; tEnd=max(tn);
    if tNov<=t0 || tEnd<=tNov, u=(tn-t0)/max(tEnd-t0,1);
    else
        u=zeros(size(tn)); pre=tn<tNov;
        u(pre)=0.27*(tn(pre)-t0)/(tNov-t0); u(~pre)=0.27+0.73*(tn(~pre)-tNov)/(tEnd-tNov);
    end
    grid=interp1(u,tv,linspace(0,1,nGrid),'linear',NaN);
end

function z = zscoreSafe(y)
    s=std(y,'omitnan'); if s<1e-9, z=y-mean(y,'omitnan'); else, z=(y-mean(y,'omitnan'))/s; end
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

function cmap = viridisMap()
    a=[68 1 84;59 82 139;33 145 140;94 201 98;253 231 37]/255;
    x=linspace(0,1,size(a,1)); xi=linspace(0,1,256);
    cmap=[interp1(x,a(:,1),xi)' interp1(x,a(:,2),xi)' interp1(x,a(:,3),xi)'];
end