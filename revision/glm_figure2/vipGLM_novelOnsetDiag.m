function vipGLM_novelOnsetDiag()
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_NOVELONSETDIAG  Diagnose the novel-onset upregulation calculation
% WITHOUT re-running the cohort GLM. Reads the saved cellTable (for each cell's
% GLM primary feature) + the session files (for per-trial dF/F), and prints,
% per cell, the last-10-familiar vs first-10-novel window means, delta, and
% one-sided paired-t p-value. Also prints the same with an UNPAIRED rank-sum so
% you can compare. Focuses on (but is not limited to) the Licking group.
%
% STEP 1: saved cohort results .mat (res.cellTable)
% STEP 2: the session .mat files (multi-select)

    uiwait(msgbox({'STEP 1 of 2:','','Select saved cohort results (res.cellTable).'},'Step 1','help','modal'));
        uiwait(msgbox({'Select the COHORT RESULTS file.','','ONE file: the cohortGLM_*.mat saved by vipGLM_cohort (has res.cellTable).','','NOT a session/m8xxx file.'},'Cohort results (ONE file)','help'));
    [rf,rp]=uigetfile('*.mat','COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes); if isequal(rf,0), return; end
    R=load(fullfile(rp,rf));
    if ~isfield(R,'res')||~isfield(R.res,'cellTable'), error('No res.cellTable.'); end
    T=R.res.cellTable; T=T(T.context=="Familiar" & T.ok,:);

    uiwait(msgbox({'STEP 2 of 2:','','Select the SESSION .mat files (multi-select).'},'Step 2','help','modal'));
        uiwait(msgbox({'Select the SESSION DATA files.','','ONE OR MORE m8xxx session .mat files (hold Cmd/Shift for several).'},'Session data (one or more)','help'));
    [sf,sp]=uigetfile('*.mat','SESSION DATA files (ONE OR MORE m8xxx .mat; multi-select)',dfSess,'MultiSelect','on'); if isequal(sf,0), return; end
    if ischar(sf), sf={sf}; end

    catNames={'Vel','Vel','Accel','Licking','Reward','Position','Pupil'};
    % tally per group: [nIncreased_paired, nIncreased_unpaired, nTotal]
    grp=containers.Map('KeyType','char','ValueType','any');
    feats={'Velocity','Acceleration','Licking','Reward','Position','Pupil'};
    for k=1:numel(feats), grp(feats{k})=[0 0 0]; end

    fprintf('\n================ NOVEL-ONSET DIAGNOSTIC ================\n');
    for si=1:numel(sf)
        fn=sf{si}; S=load(fullfile(sp,fn));
        if ~isfield(S,'sData'), warning('skip %s',fn); continue; end
        sData=S.sData;
        if ~(isfield(sData.imdata,'roiSignals')&&numel(sData.imdata.roiSignals)>=2), continue; end
        dff=sData.imdata.roiSignals(2).dff; nCells=size(dff,1); nF=size(dff,2);
        ti=double(sData.behavior.signals.trialIndexesDs); ti=ti(:); N=min(nF,numel(ti)); ti=ti(1:N);
        ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2);
        ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
        fam_id=ctx(1);
        vtr=ti>0 & ti<=length(ctx); cctx=nan(N,1); cctx(vtr)=ctx(ti(vtr)); cctx=fillmissing(cctx,'nearest');
        allTrials=unique(ti(ti>0));
        isFam=arrayfun(@(t) mode(cctx(ti==t))==fam_id, allTrials);
        novTrials=allTrials(~isFam);
        if isempty(novTrials), continue; end
        firstNovTrial=novTrials(1);
        % INITIAL familiar block only = familiar trials BEFORE the first novel trial
        % (excludes the post-novel return-to-familiar block).
        initFam=allTrials(isFam & allTrials<firstNovTrial);
        if numel(initFam)<3 || numel(novTrials)<3, continue; end
        lastFam=initFam(max(1,end-9):end); firstNov=novTrials(1:min(10,end));

        rows=T(strcmpi(string(T.session),string(fn)),:);
        for r=1:height(rows)
            cid=rows.cellID(r); if cid<1||cid>nCells, continue; end
            lab=char(string(rows.primaryFeature(r)));
            if ~isKey(grp,lab), continue; end
            y=double(dff(cid,1:N))';
            fM=arrayfun(@(t) mean(y(ti==t),'omitnan'), lastFam);
            nM=arrayfun(@(t) mean(y(ti==t),'omitnan'), firstNov);
            fM=fM(~isnan(fM)); nM=nM(~isnan(nM));
            if numel(fM)<3||numel(nM)<3, continue; end
            np=min(numel(fM),numel(nM));
            pPair=NaN; pUnp=NaN;
            try [~,pPair]=ttest(nM(1:np),fM(1:np),'Tail','right'); catch, end
            try pUnp=ranksum(nM,fM,'tail','right'); catch, end
            incP=pPair<0.05; incU=pUnp<0.05;
            v=grp(lab); v=v+[incP incU 1]; grp(lab)=v;
            if strcmp(lab,'Licking')
                fprintf('LICK %s c%-3d famTrials[%d..%d] novTrials[%d..%d]\n', ...
                    fn(1:min(8,end)), cid, lastFam(1), lastFam(end), firstNov(1), firstNov(end));
                fprintf('     famVals=%s\n', mat2str(round(fM(:)',3)));
                fprintf('     novVals=%s\n', mat2str(round(nM(:)',3)));
                fprintf('     famM=%+.3f novM=%+.3f d=%+.3f pPair=%.3f pUnp=%.3f\n', ...
                    mean(fM), mean(nM), mean(nM)-mean(fM), pPair, pUnp);
            end
        end
    end

    fprintf('\n--------- fraction increased per group (paired vs unpaired) ---------\n');
    for k=1:numel(feats)
        v=grp(feats{k});
        if v(3)>0
            fprintf('%-13s  paired %3.0f%%   unpaired %3.0f%%   (n=%d)\n', ...
                feats{k}, 100*v(1)/v(3), 100*v(2)/v(3), v(3));
        end
    end
    % overall
    tot=[0 0 0]; ks=keys(grp); for i=1:numel(ks), tot=tot+grp(ks{i}); end
    if tot(3)>0
        fprintf('%-13s  paired %3.0f%%   unpaired %3.0f%%   (n=%d)\n', ...
            'ALL', 100*tot(1)/tot(3), 100*tot(2)/tot(3), tot(3));
    end
    fprintf('=======================================================\n');
end
