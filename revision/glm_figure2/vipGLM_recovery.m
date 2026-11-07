function vipGLM_recovery()
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_RECOVERY  Ground-truth validation of the encoding-model pipeline.
%
% AIM: prove the GLM recovers known tuning. Builds surrogate cells with KNOWN
% injected tuning from the REAL behavior of a chosen session (preserving the
% true collinearity between velocity / licking / position / reward / pupil),
% then checks whether the GLM correctly recovers the injected feature set.
% This is the methods-validation figure (shows the pipeline isn't inventing
% structure). Self-contained: fitting math matches vipGLM_fitCell.
%
% INPUT:  ONE SESSION DATA file (sData .mat) to draw real behavior from, via prompt.
% OUTPUT: recovery-performance figure(s). Requires vipGLM_defaultFolder.m.
%
%   Self-contained: no dependency on vipGLM_fitCell.m. Fitting math is identical
%   (causal shifts, RBF position, context z-scoring, manual ridge with
%   unpenalized intercept, pooled-fold CV R^2, raw-variable circular-shift
%   null) but the design matrix is built ONCE per surrogate and only the
%   shuffled block is rebuilt inside the permutation loop. Conditions run
%   under parfor.
%
%   HONEST SCOPE: a Poisson-spikes-through-softplus generator validates the
%   GLM as an approximate LINEAR readout of WHICH variables drive a cell, not
%   as recovery of the true nonlinear encoding function.

    userHome = getenv('HOME');
    folderPath = fullfile(userHome,'Desktop','VIP-GCaMP-sData');
    if ~exist(folderPath,'dir'), folderPath = fullfile(userHome,'Desktop'); end
        uiwait(msgbox({'Select ONE SESSION DATA file.','','A single m8xxx session .mat file.'},'Session data (one file)','help'));
    [fileName,pathName] = uigetfile('*.mat','ONE SESSION DATA file (m8xxx .mat) for REAL behavior',dfSess);
    if isequal(fileName,0), return; end

    uiFig = uifigure('Name',['GLM Recovery Validator (fast): ' fileName],'Position',[100 100 380 490]);
    g = uigridlayout(uiFig,[15 2]); g.RowHeight=repmat({'fit'},1,15); g.ColumnWidth={'fit','1x'};

    uilabel(g,'Text','Context:');                 envDD = uidropdown(g,'Items',{'Familiar Context','Novel Context'},'Value','Novel Context');
    uilabel(g,'Text','Significance:');             sigDD = uidropdown(g,'Items',{'Permutation (rigorous)','t-test (fast)'},'Value','Permutation (rigorous)');
    uilabel(g,'Text','Lambda mode:');              lamDD = uidropdown(g,'Items',{'Fixed lambda','Nested CV lambda'},'Value','Fixed lambda');
    uilabel(g,'Text','Ridge \lambda:');           lamF  = uieditfield(g,'numeric','Value',0.5);
    uilabel(g,'Text','# permutations:');          permF = uieditfield(g,'numeric','Value',200,'Limits',[20 2000],'RoundFractionalValues','on');
    uilabel(g,'Text','Lag min (s):');             lagMinF = uieditfield(g,'numeric','Value',-3,'Limits',[-10 10]);
    uilabel(g,'Text','Lag max (s):');             lagMaxF = uieditfield(g,'numeric','Value',3,'Limits',[-10 10]);
    uilabel(g,'Text','# lag points:');            nLagF   = uieditfield(g,'numeric','Value',5,'Limits',[1 21],'RoundFractionalValues','on');
    uilabel(g,'Text','# example cells:');         nExF  = uieditfield(g,'numeric','Value',8,'Limits',[2 50],'RoundFractionalValues','on');
    uilabel(g,'Text','SNR sweep points:');        snrF  = uieditfield(g,'numeric','Value',5,'Limits',[1 10],'RoundFractionalValues','on');
    uilabel(g,'Text','Repeats per condition:');   repF  = uieditfield(g,'numeric','Value',10,'Limits',[1 50],'RoundFractionalValues','on');
    uilabel(g,'Text','Random seed:');             seedF = uieditfield(g,'numeric','Value',42);
    uilabel(g,'Text','Use parfor:');              parDD = uidropdown(g,'Items',{'Yes','No'},'Value','Yes');
    uilabel(g,'Text','Save results .mat:');       saveDD= uidropdown(g,'Items',{'Yes','No'},'Value','Yes');

    runB = uibutton(g,'Text','Run recovery validation', 'ButtonPushedFcn', @(s,e) runValidation( ...
        fullfile(pathName,fileName), envDD.Value, lamF.Value, round(permF.Value), round(nExF.Value), ...
        round(snrF.Value), round(repF.Value), seedF.Value, strcmp(parDD.Value,'Yes'), strcmp(saveDD.Value,'Yes'), ...
        contains(sigDD.Value,'t-test'), contains(lamDD.Value,'Nested'), lagMinF.Value, lagMaxF.Value, round(nLagF.Value), uiFig));
    runB.Layout.Column=[1 2];
end

% =====================================================================
function runValidation(filePath, envChoice, lambdaVal, nPerms, nExample, nSNR, nRep, seed, usePar, doSave, useTtest, useNested, lagMin, lagMax, nLags, uiFig)
    rng(seed);
    d = uiprogressdlg(uiFig,'Title','Recovery validation','Message','Loading session...','Value',0.02);
    cleanup = onCleanup(@() closeIfValid(d)); %#ok<NASGU>

    D = load(filePath);
    if ~isfield(D,'sData'), uialert(uiFig,'No sData.','Error'); return; end
    sData = D.sData;
    if ~(isfield(sData.imdata,'roiSignals') && numel(sData.imdata.roiSignals)>=2)
        uialert(uiFig,'No roiSignals(2).dff.','Error'); return; end
    dff_all = sData.imdata.roiSignals(2).dff;
    nRealCells = size(dff_all,1);

    d.Message='Extracting behavior & building design scaffold...'; d.Value=0.06;
    B = extractBehavior(sData, size(dff_all,2));
    Fs=B.Fs; N=B.N;
    mask = contextMask(sData,B,envChoice,N);
    if sum(mask)<100, uialert(uiFig,'Too few frames in context.','Error'); return; end

    S = buildScaffold(B, mask, Fs, lagMin, lagMax, nLags);

    d.Message='Estimating R^2 / reliability from real cells...'; d.Value=0.10;
    exCells = unique(round(linspace(1,nRealCells,min(nExample,nRealCells))));
    realR2=nan(numel(exCells),1); realRel=nan(numel(exCells),1); realNoise=nan(numel(exCells),1);
    for k=1:numel(exCells)
        raw = double(dff_all(exCells(k),1:N))';
        rr = raw(mask);                       % mask first
        ym = S.ctx_z(rr); ym(isnan(ym))=0;    % then within-context z-score
        realR2(k) = pooledRidgeR2(S.Xfull_m, ym, S.cv_idx, S.k_folds, lambdaVal);
        realRel(k)= reliabilityStat(raw, B, mask);
        % dF/F noise floor: high-frequency component of the RAW (un-z-scored) trace
        realNoise(k) = std(diff(rr),'omitnan')/sqrt(2);   % first-difference noise estimate
    end
    realR2=realR2(~isnan(realR2));
    if isempty(realR2), target_R2=0.2; R2_band=[0.1 0.3];
    else, target_R2=median(realR2); R2_band=prctile(realR2,[25 75]); end
    relTarget = median(realRel,'omitnan'); if isnan(relTarget), relTarget=0.6; end
    noiseFloor = median(realNoise,'omitnan');   % anchors surrogate measurement noise

    feat=S.feat; nf=numel(feat); rawF={B.vel,B.acc,B.lick,B.valve,B.pos,B.pupil};
    truth={}; labels={};
    for fI=1:nf, t=false(1,nf); t(fI)=true; truth{end+1}=t; labels{end+1}=['Pure: ' feat{fI}]; end %#ok<AGROW>
    truth{end+1}=ismember(feat,{'Velocity','Licking'});         labels{end+1}='Combo: Vel+Lick';
    truth{end+1}=ismember(feat,{'Position','Reward'});           labels{end+1}='Combo: Pos+Reward';
    truth{end+1}=ismember(feat,{'Velocity','Position','Pupil'}); labels{end+1}='Combo: Vel+Pos+Pupil';
    truth{end+1}=ismember(feat,{'Velocity'});                    labels{end+1}='TRAP: Vel only';
    truth{end+1}=false(1,nf);                                    labels{end+1}='NOISE: none';
    if ~B.hasPupil
        keep=true(1,numel(truth)); for c=1:numel(truth), if truth{c}(6), keep(c)=false; end; end
        truth=truth(keep); labels=labels(keep);
    end
    nCond=numel(truth);

    alphas = linspace(0.3, 2.5, nSNR);   % modulation-strength sweep (difficulty axis)

    d.Message='Fitting surrogates...'; d.Value=0.15;
    sigTensor=false(nCond,nf,nSNR,nRep); r2Tensor=nan(nCond,nSNR,nRep);

    [Ci,Gi,Ri]=ndgrid(1:nCond,1:nSNR,1:nRep);
    Ci=Ci(:); Gi=Gi(:); Ri=Ri(:); nJobs=numel(Ci);
    sigFlat=false(nJobs,nf); r2Flat=nan(nJobs,1);

    truthLoc=truth; rawFLoc=rawF; SLoc=S; BLoc=B; maskLoc=mask; lam=lambdaVal; nP=nPerms; rel=relTarget; nfloor=noiseFloor; tt=useTtest; nn=useNested;
    seeds = seed + (1:nJobs);

    runPar = usePar && ~isempty(ver('parallel'));
    if runPar
        % parfor body must be transparent: no nested functions, no UI access
        % inside the loop. Progress is reported via a DataQueue + a handle-class
        % counter (progressTick, defined at file end) updated on the CLIENT.
        dq = parallel.pool.DataQueue;
        progressTick(d, nJobs, true);   % reset counter
        afterEach(dq, @(~) progressTick(d, nJobs, false));
        d.Message = sprintf('Fitting %d surrogates on pool...', nJobs);
        parfor j=1:nJobs
            rng(seeds(j));
            w=zeros(1,nf); w(truthLoc{Ci(j)})=1;
            ys = makeSurrogateScaffold(rawFLoc, w, alphas(Gi(j)), BLoc, maskLoc, SLoc, rel, nfloor);
            [r2j,sigj] = fitLean(ys, SLoc, maskLoc, rawFLoc, BLoc, lam, nP, tt, nn);
            sigFlat(j,:)=sigj; r2Flat(j)=r2j;
            send(dq,j);
        end
    else
        for j=1:nJobs
            rng(seeds(j));
            w=zeros(1,nf); w(truthLoc{Ci(j)})=1;
            ys = makeSurrogateScaffold(rawFLoc, w, alphas(Gi(j)), BLoc, maskLoc, SLoc, rel, nfloor);
            [r2j,sigj] = fitLean(ys, SLoc, maskLoc, rawFLoc, BLoc, lam, nP, tt, nn);
            sigFlat(j,:)=sigj; r2Flat(j)=r2j;
            if mod(j,5)==0||j==nJobs, d.Value=0.15+0.8*(j/nJobs); d.Message=sprintf('Fitting %d/%d...',j,nJobs); end
        end
    end

    for j=1:nJobs
        sigTensor(Ci(j),:,Gi(j),Ri(j))=sigFlat(j,:);
        r2Tensor(Ci(j),Gi(j),Ri(j))=r2Flat(j);
    end

    medR2byGain=squeeze(median(median(r2Tensor,3,'omitnan'),1,'omitnan'));
    [~,giMatch]=min(abs(medR2byGain-target_R2));

    d.Value=0.97; d.Message='Plotting...';
    plotRecovery(sigTensor,r2Tensor,truth,labels,feat,alphas,giMatch,R2_band,target_R2,medR2byGain,B,mask,envChoice,lambdaVal,nPerms,useTtest,useNested);

    % Seeing is believing: show example surrogate dF/F traces (at matched alpha)
    % next to a real example cell, in the style of vipGLM_inspect.
    plotSurrogateExamples(rawF, truth, labels, alphas(giMatch), B, mask, S, relTarget, noiseFloor, ...
        dff_all, exCells, envChoice);

    if doSave
        outFile=sprintf('recovery_%s_%s.mat', erase(envChoice,' '), datestr(now,'yyyymmdd_HHMMSS'));
        res=struct('sigTensor',sigTensor,'r2Tensor',r2Tensor,'truth',{truth},'labels',{labels}, ...
            'feat',{feat},'alphas',alphas,'giMatch',giMatch,'R2_band',R2_band,'target_R2',target_R2, ...
            'medR2byGain',medR2byGain,'meta',struct('env',envChoice,'lambda',lambdaVal, ...
            'lambdaMode',ternaryStr(useNested,'nested','fixed'),'sigTest',ternaryStr(useTtest,'ttest','perm'), ...
            'nPerms',nPerms,'lagMin',lagMin,'lagMax',lagMax,'nLags',nLags, ...
            'nRep',nRep,'seed',seed,'file',filePath)); %#ok<NASGU>
        save(outFile,'res'); fprintf('Saved %s\n',outFile);
    end

    fprintf('\n=== Recovery validation ===\nReal-cell R^2 median %.2f IQR [%.2f %.2f] (n=%d)\n', ...
        target_R2,R2_band(1),R2_band(2),numel(realR2));
    fprintf('Matched slice: alpha=%.2f surrogate median R^2=%.2f\n',alphas(giMatch),medR2byGain(giMatch));
end

% =====================================================================
function S = buildScaffold(B, mask, Fs, lagMin, lagMax, nLags)
    S.feat={'Velocity','Acceleration','Licking','Reward','Position','Pupil'};
    % within-context z-score: mask is the context, so operate on masked vectors directly
    zc = @(v) (v-mean(v,'omitnan'))./max(std(v,'omitnan'),1e-8);
    S.ctx_z = zc;   % stored for the generator (applied to masked vectors)
    sigma_frames=round(0.125*Fs); S.kernel=gausswin(round(sigma_frames*3)); S.kernel=S.kernel/sum(S.kernel);
    S.shifts=round(linspace(lagMin,lagMax,nLags)*Fs);   % lag window [lagMin,lagMax], nLags points
    S.num_bases=15; S.pos_centers=linspace(0,max(B.pos),S.num_bases); S.pos_width=(max(B.pos)/S.num_bases)*1.2;

    % masked (within-context) raw behaviors
    S.vel_m=B.vel(mask); S.acc_m=B.acc(mask); S.lick_m=B.lick(mask);
    S.valve_m=B.valve(mask); S.pos_m=B.pos(mask); S.pupil_m=B.pupil(mask);
    S.raw_masked={S.vel_m,S.acc_m,S.lick_m,S.valve_m,S.pos_m,S.pupil_m};
    Nm=sum(mask); S.Nm=Nm;

    Xv =mkShift(conv(zc(S.vel_m),  S.kernel,'same'),S.shifts);
    Xa =mkShift(conv(zc(S.acc_m),  S.kernel,'same'),S.shifts);
    Xl =mkShift(conv(zc(S.lick_m), S.kernel,'same'),S.shifts);
    Xr =mkShift(conv(zc(S.valve_m),S.kernel,'same'),S.shifts);
    Xp =zeros(Nm,S.num_bases); for b=1:S.num_bases, Xp(:,b)=exp(-((S.pos_m-S.pos_centers(b)).^2)/(2*S.pos_width^2)); end
    Xpu=mkShift(conv(zc(S.pupil_m),S.kernel,'same'),S.shifts);

    ns=numel(S.shifts);
    S.idx_vel=1:ns; S.idx_acc=ns+(1:ns); S.idx_lick=2*ns+(1:ns); S.idx_rew=3*ns+(1:ns);
    S.idx_pos=4*ns+(1:S.num_bases); S.idx_pupil=4*ns+S.num_bases+(1:ns);
    S.feat_indices={S.idx_vel,S.idx_acc,S.idx_lick,S.idx_rew,S.idx_pos,S.idx_pupil};

    S.Xfull_m=[Xv,Xa,Xl,Xr,Xp,Xpu]; S.Xfull_m(isnan(S.Xfull_m))=0;

    S.k_folds=10; cs=round(Fs*10); nc=ceil(Nm/cs);
    fa=repmat(1:S.k_folds,1,ceil(nc/S.k_folds)); cv=zeros(Nm,1);
    for i=1:nc, cv(((i-1)*cs+1):min(i*cs,Nm))=fa(i); end
    S.cv_idx=cv;
end

function [R2_full, sig] = fitLean(ym, S, mask, rawF, B, lambdaVal, nPerms, useTtest, useNested)
    if nargin<8, useTtest=false; end
    if nargin<9, useNested=false; end
    nf=6; sig=false(1,nf);
    ym=ym(:); ym(isnan(ym))=0;             % already within-context (masked length)
    Nm=S.Nm; zc=S.ctx_z;

    % resolve lambda: fixed, or CV-tuned grid search (matches vipGLM_fitCell)
    % (NB: tuned on the same folds used for evaluation -- not true nested CV)
    if useNested
        grid=[0.01 0.1 0.5 1 2 5 10 50]; gr=zeros(numel(grid),1);
        for gi=1:numel(grid), gr(gi)=pooledRidgeR2(S.Xfull_m,ym,S.cv_idx,S.k_folds,grid(gi)); end
        [~,bi]=max(gr); lambdaVal=grid(bi);
    end

    R2_full = pooledRidgeR2(S.Xfull_m, ym, S.cv_idx, S.k_folds, lambdaVal);
    if R2_full<=0, return; end

    dR2=zeros(1,nf); pv=ones(1,nf);
    for f=1:nf
        keep=setdiff(1:size(S.Xfull_m,2),S.feat_indices{f});
        Xred_m=S.Xfull_m(:,keep);

        if useTtest
            % Fast: paired t-test on per-fold MSE (full vs reduced), no permutations
            [R2_red, mse_red, mse_full] = foldMSEcompare(Xred_m, S.Xfull_m, ym, S.cv_idx, S.k_folds, lambdaVal);
            dR2(f)=R2_full-R2_red;
            try [~,pv(f)]=ttest(mse_red, mse_full, 'Tail','right'); catch, pv(f)=1; end
        else
            % Rigorous: within-context circular-shift permutation null
            R2_red=pooledRidgeR2(Xred_m,ym,S.cv_idx,S.k_folds,lambdaVal);
            dR2(f)=R2_full-R2_red;
            isPos=strcmp(S.feat{f},'Position'); rv=S.raw_masked{f};   % within-context raw vector
            ndr=zeros(nPerms,1);
            for p=1:nPerms
                sa=randi([round(Nm*0.1),round(Nm*0.9)]); rs=circshift(rv,sa,1);
                if ~isPos
                    blk=mkShift(conv(zc(rs),S.kernel,'same'),S.shifts);
                else
                    blk=zeros(Nm,S.num_bases);
                    for b=1:S.num_bases, blk(:,b)=exp(-((rs-S.pos_centers(b)).^2)/(2*S.pos_width^2)); end
                end
                Xs=S.Xfull_m; Xs(:,S.feat_indices{f})=blk; Xs(isnan(Xs))=0;
                ndr(p)=pooledRidgeR2(Xs,ym,S.cv_idx,S.k_folds,lambdaVal)-R2_red;
            end
            pv(f)=(sum(ndr>=dR2(f))+1)/(nPerms+1);
        end
    end
    sig=(pv<0.05)&(dR2>0.01);
    if ~B.hasPupil, sig(6)=false; end
end

% Per-fold MSE for full and reduced models (for the t-test significance mode)
function [pooled_R2_red, mse_red, mse_full] = foldMSEcompare(Xred, Xfull, y, cv_idx, k_folds, lambdaVal)
    yp_red=nan(length(y),1); mse_red=nan(k_folds,1); mse_full=nan(k_folds,1);
    for f=1:k_folds
        te=(cv_idx==f); tr=(cv_idx~=f & cv_idx>0);
        if ~any(tr)||~any(te), continue; end
        Xtr_r=[ones(sum(tr),1),Xred(tr,:)]; Xte_r=[ones(sum(te),1),Xred(te,:)];
        pr=lambdaVal*eye(size(Xtr_r,2)); pr(1,1)=0; br=(Xtr_r'*Xtr_r+pr)\(Xtr_r'*y(tr));
        yr=Xte_r*br; yp_red(te)=yr;
        Xtr_f=[ones(sum(tr),1),Xfull(tr,:)]; Xte_f=[ones(sum(te),1),Xfull(te,:)];
        pf=lambdaVal*eye(size(Xtr_f,2)); pf(1,1)=0; bf=(Xtr_f'*Xtr_f+pf)\(Xtr_f'*y(tr));
        yf=Xte_f*bf;
        mse_red(f)=mean((y(te)-yr).^2); mse_full(f)=mean((y(te)-yf).^2);
    end
    vi=~isnan(yp_red)&cv_idx>0;
    sst=sum((y(vi)-mean(y(vi))).^2); sse=sum((y(vi)-yp_red(vi)).^2);
    if sst>0, pooled_R2_red=1-sse/sst; else, pooled_R2_red=0; end
end

function R2 = pooledRidgeR2(X,y,cv_idx,k_folds,lambdaVal)
    N=length(y); yp=nan(N,1);
    for f=1:k_folds
        te=(cv_idx==f); tr=(cv_idx~=f & cv_idx>0);
        if ~any(tr)||~any(te), continue; end
        Xtr=[ones(sum(tr),1),X(tr,:)]; Xte=[ones(sum(te),1),X(te,:)];
        pen=lambdaVal*eye(size(Xtr,2)); pen(1,1)=0;
        b=(Xtr'*Xtr+pen)\(Xtr'*y(tr)); yp(te)=Xte*b;
    end
    vi=~isnan(yp)&cv_idx>0;
    if ~any(vi), R2=0; return; end
    sst=sum((y(vi)-mean(y(vi))).^2); sse=sum((y(vi)-yp(vi)).^2);
    if sst>0, R2=1-sse/sst; else, R2=0; end
end

% =====================================================================
function ys = makeSurrogateScaffold(rawF, w, alpha, B, mask, S, relTarget, noiseFloor)
% Generative surrogate for a tonically-active inhibitory neuron:
%   - high background firing rate, modulated MULTIPLICATIVELY (up & down) by
%     a feature drive: rate = baseRate * exp(alpha * z(drive))
%   - small per-spike dF/F amplitude (A_spike), so single spikes are invisible
%     and the signal emerges only through summation of the high background rate
%   - bi-exponential GCaMP6s-like kernel, normalized to UNIT PEAK so each spike
%     contributes exactly A_spike at its peak
%   - additive measurement noise set from the real-cell dF/F noise floor
%   'alpha' is the difficulty axis (modulation strength), swept by the caller.
    baseRate = 6;        % Hz, tonic background (VIP/IS-cell regime)
    A_spike  = 0.02;     % peak dF/F per single spike (measured)
    tau_rise = 0.07;     % s
    tau_dec  = 1.0;      % s, GCaMP6s decay

    Nm=S.Nm; zc=S.ctx_z; drive=zeros(Nm,1);
    rawM=S.raw_masked;
    for k=1:numel(w), if w(k)~=0, drive=drive+w(k)*zc(rawM{k}); end; end
    % normalize the combined drive to unit SD so alpha has a consistent meaning
    sd=std(drive,'omitnan'); if sd>0, drive=drive/sd; end

    % reliability knob: slow trial-to-trial jitter added to the drive
    if isempty(relTarget)||isnan(relTarget), relTarget=0.6; end
    jitterSD=max(0,(1-relTarget));
    drive=drive+jitterSD*smoothdata(randn(Nm,1),'gaussian',round(B.Fs*1.0));

    % multiplicative, bidirectional rate modulation around the baseline
    rate = baseRate .* exp(alpha .* drive);
    rate(~isfinite(rate)) = baseRate;
    rate = min(rate, 20*baseRate);          % guard against runaway exp tails

    dt=1/B.Fs; spikes=poissrnd(rate*dt);

    % unit-PEAK bi-exponential kernel, scaled by per-spike amplitude
    tk=(0:round(5*tau_dec*B.Fs))'/B.Fs;
    ck=(1-exp(-tk/tau_rise)).*exp(-tk/tau_dec);
    ck=ck/max(ck);                          % unit peak (NOT unit area)
    ca=conv(spikes,ck,'full'); ca=ca(1:Nm);
    ca=A_spike*ca;

    % measurement noise from the real-cell dF/F noise floor (physically anchored)
    if isempty(noiseFloor)||isnan(noiseFloor)||noiseFloor<=0
        nf=0.5*A_spike;                     % fallback if no estimate available
    else
        nf=noiseFloor;
    end
    ys = ca + nf*randn(Nm,1);
    % baseline-subtract (F0 as a low percentile) so surrogate dF/F sits near zero
    % like real baseline-corrected traces; harmless for the GLM (z-scored downstream)
    ys = ys - prctile(ys,8);
end

% =====================================================================
function X=mkShift(vec,shifts)
    N=length(vec); X=zeros(N,length(shifts));
    for i=1:length(shifts)
        s=shifts(i);
        if s==0, X(:,i)=vec; elseif s>0, X(:,i)=[zeros(s,1);vec(1:end-s)];
        else, sa=abs(s); X(:,i)=[vec(sa+1:end);zeros(sa,1)]; end
    end
end

function rel = reliabilityStat(dff,B,mask)
    rel=NaN;
    try
        pos=B.pos; tr=B.trialIdx; utr=unique(tr(mask&tr>0));
        if numel(utr)<4, return; end
        odd=utr(1:2:end); even=utr(2:2:end);
        nb=50; edges=linspace(0,max(pos),nb+1);
        mo=bm(dff,pos,tr,odd,edges,mask); me=bm(dff,pos,tr,even,edges,mask);
        gd=~isnan(mo)&~isnan(me); if nnz(gd)>5, rel=corr(mo(gd),me(gd)); end
    catch, end
end
function m=bm(dff,pos,tr,trials,edges,mask)
    nb=numel(edges)-1; m=nan(nb,1); sel=mask&ismember(tr,trials);
    for b=1:nb, idx=sel&pos>=edges(b)&pos<edges(b+1); if any(idx), m(b)=mean(dff(idx),'omitnan'); end; end
end

% =====================================================================
function B=extractBehavior(sData,dffLen)
    time_len=length(sData.behavior.signals.timeDs);
    pos_len=length(sData.behavior.signals.corridorPositionDs);
    vel_len=length(sData.behavior.signals.velocityDs);
    N=min([dffLen,time_len,pos_len,vel_len]);
    time=sData.behavior.signals.timeDs(1:N); time=time(:); Fs=1/mean(diff(time),'omitnan');
    pos=double(sData.behavior.signals.corridorPositionDs(1:N)); pos=pos(:);
    vel=double(sData.behavior.signals.velocityDs(1:N)); vel=vel(:);
    lick=double(sData.behavior.signals.lickIfreqDs(1:N)); lick=lick(:);
    if length(vel)>2, vel(1:2)=vel(3); end
    acc=gradient(vel)*Fs; acc=acc(:);
    valve=zeros(N,1);
    if isfield(sData.behavior,'signals')&&isfield(sData.behavior.signals,'rewardDs')
        vt=double(sData.behavior.signals.rewardDs); valve(1:min(N,length(vt)))=vt(1:min(N,length(vt)))>0;
    elseif isfield(sData,'daqdata')&&isfield(sData.daqdata,'frameIndex')&&isfield(sData.daqdata,'waterValve')
        vr=double(sData.daqdata.waterValve>0.5); ft=find(diff(sData.daqdata.frameIndex)==1)+1; lim=min(N,length(ft)-1);
        for i=1:lim, if any(vr(ft(i):ft(i+1)-1)), valve(i)=1; end; end
        if N<=length(ft)&&any(vr(ft(N):end)), valve(N)=1; end
    end
    pupil=nan(N,1); hasPupil=false;
    try
        if isfield(sData.behavior,'pupildata')&&isfield(sData.behavior.pupildata,'diameter')
            pd=double(sData.behavior.pupildata.diameter); pd(isnan(pd))=inf;
            pfi=find(diff(sData.behavior.signals.pupilCamFrameIndexCorr)==1)+1; ifi=find(diff(sData.daqdata.frameIndex)==1);
            pUs=nan(size(sData.behavior.signals.pupilCamFrameIndexCorr)); vl=min(length(pfi),length(pd)); pUs(pfi(1:vl))=pd(1:vl);
            pUs=fillmissing(pUs,'nearest'); pUs(isinf(pUs))=nan; pr=double(pUs(ifi));
            if length(pr)>=N, pupil=pr(1:N); else, pupil=[pr(:);nan(N-length(pr),1)]; end
            pupil=pupil(:); hasPupil=~all(isnan(pupil));
        end
    catch, end
    if hasPupil, pupil=fillmissing(pupil,'nearest'); pupil(isnan(pupil))=mean(pupil,'omitnan'); else, pupil=zeros(N,1); end
    trialIdx=sData.behavior.signals.trialIndexesDs(1:N); trialIdx=trialIdx(:);
    B=struct('N',N,'Fs',Fs,'pos',pos,'vel',vel,'acc',acc,'lick',lick,'valve',valve,'pupil',pupil,'hasPupil',hasPupil,'trialIdx',trialIdx);
end

function mask=contextMask(sData,B,envChoice,N)
    ctx_raw=mode(sData.behavior.trialMatrices.contextMatrix,2); ctx=round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
    tr=B.trialIdx; vtr=tr>0&tr<=length(ctx); cctx=nan(N,1); cctx(vtr)=ctx(tr(vtr)); cctx=fillmissing(cctx,'nearest');
    fam=ctx(1); u=unique(ctx); nov=u(u~=fam);
    if strcmp(envChoice,'Familiar Context')
        mask=(cctx==fam);   % all familiar frames (both blocks)
    else
        if isempty(nov), mask=false(N,1); else, mask=(cctx==nov(1)); end
    end
    mask=mask(:);
end

% =====================================================================
function plotSurrogateExamples(rawF, truth, labels, alpha, B, mask, S, relTarget, noiseFloor, dff_all, exCells, envChoice)
% Show example surrogate dF/F traces (trace + spatial average) next to one
% real cell, so the surrogates can be eyeballed against real VIP data.
    nf=numel(S.feat);
    % choose a few representative conditions to display
    wantLabels = {'Pure: Velocity','Pure: Position','Combo: Vel+Lick','NOISE: none'};
    showIdx = [];
    for i=1:numel(wantLabels)
        h=find(strcmp(labels,wantLabels{i}),1);
        if ~isempty(h), showIdx(end+1)=h; end %#ok<AGROW>
    end
    nShow = numel(showIdx);
    Fs=B.Fs; N=B.N; pos=B.pos;
    pos_m = pos(mask);                        % within-context position

    % a real example cell for comparison (raw dF/F), masked to the context
    realCell = exCells(max(1,round(numel(exCells)/2)));
    realDff_m = double(dff_all(realCell,1:N))'; realDff_m = realDff_m(mask);

    % spatial-average helper operates on masked-length vectors
    numBins=50; edges=linspace(0,max(pos),numBins+1); ctrs=linspace(0,max(pos),numBins);
    binIt=@(v) arrayfun(@(b) mean(v(pos_m>=edges(b) & pos_m<edges(b+1)),'omitnan'), 1:numBins);

    f=figure('Color','w','Position',[60 60 1400 900],'Name','Surrogate dF/F examples');
    nRows=nShow+1;           % +1 for the real cell
    tWin=60;                 % seconds shown in the trace column
    tvec=(1:S.Nm)'/Fs;       % masked-context time axis (starts at 0 within the context)
    win = tvec<=tWin;        % first 60 s of the context

    for r=1:nRows
        if r<=nShow
            w=zeros(1,nf); w(truth{showIdx(r)})=1;
            ys = makeSurrogateScaffold(rawF, w, alpha, B, mask, S, relTarget, noiseFloor);
            rowLabel = labels{showIdx(r)};
            col=[0.15 0.35 0.6];
        else
            ys = realDff_m;
            rowLabel = sprintf('REAL cell #%d', realCell);
            col=[0 0 0];
        end

        % trace column
        ax=subplot(nRows,3,(r-1)*3+[1 2]); hold on;
        plot(tvec(win), ys(win), 'Color', col, 'LineWidth', 0.8);
        ylabel('\DeltaF/F'); title(rowLabel,'Interpreter','none');
        if r==nRows, xlabel('Time (s)'); else, set(ax,'XTickLabel',[]); end
        axis tight; set(ax,'TickDir','out'); box off;

        % spatial-average column
        ax2=subplot(nRows,3,(r-1)*3+3); hold on;
        m=binIt(ys);
        plot(ctrs, m, 'Color', col, 'LineWidth', 1.5);
        xline([100 220],'--r','Alpha',0.4);
        if r==nRows, xlabel('Position (cm)'); else, set(ax2,'XTickLabel',[]); end
        ylabel('mean \DeltaF/F'); set(ax2,'TickDir','out'); box off; axis tight;
    end
    sgtitle(sprintf('Surrogate \\DeltaF/F examples vs real cell | %s | \\alpha=%.2f | baseline 6 Hz, 0.02 \\DeltaF/F per spike', ...
        envChoice, alpha),'FontWeight','bold');
end


% =====================================================================
function plotRecovery(sigTensor,r2Tensor,truth,labels,feat,gains,giMatch,R2_band,target_R2,medR2byGain,B,mask,envChoice,lambdaVal,nPerms,useTtest,useNested)
    nCond=size(sigTensor,1); nf=numel(feat); nSNR=size(sigTensor,3);
    TM=cell2mat(truth(:));
    detMatch=squeeze(mean(sigTensor(:,:,giMatch,:),4));

    figure('Color','w','Position',[40 40 1500 900],'Name','GLM recovery validation');

    ax1=subplot(2,3,1); imagesc(detMatch,[0 1]); colormap(ax1,parula); cb=colorbar; cb.Label.String='P(flagged)';
    set(ax1,'XTick',1:nf,'XTickLabel',feat,'XTickLabelRotation',25,'YTick',1:nCond,'YTickLabel',labels,'TickDir','out');
    title(sprintf('Recovery @ matched SNR (R^2\\approx%.2f)',medR2byGain(giMatch))); hold on;
    for c=1:nCond, for k=1:nf, if TM(c,k), rectangle('Position',[k-0.5 c-0.5 1 1],'EdgeColor',[0 0.8 0],'LineWidth',2); end; end; end
    xlabel('GLM-flagged feature'); ylabel('Surrogate (green = injected)');

    subplot(2,3,2); hold on; hits=detMatch(TM); fas=detMatch(~TM);
    histogram(hits,0:0.1:1,'FaceColor',[0.2 0.6 0.3],'FaceAlpha',0.7);
    histogram(fas,0:0.1:1,'FaceColor',[0.8 0.3 0.3],'FaceAlpha',0.6);
    xlabel('P(flagged) @ matched SNR'); ylabel('# feature-cells');
    legend({'Injected (hit)','Not injected (FA)'},'Location','north');
    title(sprintf('Hit %.2f  vs  FA %.2f',mean(hits,'omitnan'),mean(fas,'omitnan'))); set(gca,'TickDir','out'); box off;

    cTrap=find(strcmp(labels,'TRAP: Vel only'));
    subplot(2,3,3);
    if ~isempty(cTrap)
        det=detMatch(cTrap,:); b=bar(1:nf,det); b.FaceColor='flat';
        for k=1:nf, if TM(cTrap,k), b.CData(k,:)=[0.2 0.6 0.3]; else, b.CData(k,:)=[0.8 0.4 0.4]; end; end
        set(gca,'XTick',1:nf,'XTickLabel',feat,'XTickLabelRotation',25,'TickDir','out'); ylim([0 1]); ylabel('P(flagged)');
        title('Collinearity trap: Velocity only'); yline(0.05,'k--','chance');
    else, text(0.5,0.5,'No trap condition','HorizontalAlignment','center'); end
    box off;

    subplot(2,3,4); hold on; hitByG=nan(1,nSNR); faByG=nan(1,nSNR);
    for gi=1:nSNR, dm=squeeze(mean(sigTensor(:,:,gi,:),4)); hitByG(gi)=mean(dm(TM),'omitnan'); faByG(gi)=mean(dm(~TM),'omitnan'); end
    plot(medR2byGain,hitByG,'-o','Color',[0.2 0.6 0.3],'LineWidth',1.5);
    plot(medR2byGain,faByG,'-o','Color',[0.8 0.3 0.3],'LineWidth',1.5);
    yl=[0 1]; patch([R2_band(1) R2_band(2) R2_band(2) R2_band(1)],[yl(1) yl(1) yl(2) yl(2)],[0.6 0.6 0.6],'FaceAlpha',0.15,'EdgeColor','none');
    xline(target_R2,'k--'); xlabel('Surrogate R^2'); ylabel('Rate'); ylim(yl);
    legend({'Hit','FA','Real-cell R^2 band'},'Location','east'); title('Recovery vs SNR'); set(gca,'TickDir','out'); box off;

    subplot(2,3,5);
    M=[B.vel(mask) B.acc(mask) B.lick(mask) B.valve(mask) B.pos(mask) B.pupil(mask)];
    C=corr(M,'rows','pairwise'); imagesc(C,[-1 1]); colormap(gca,redbluemap()); colorbar;
    set(gca,'XTick',1:nf,'XTickLabel',feat,'XTickLabelRotation',25,'YTick',1:nf,'YTickLabel',feat,'TickDir','out');
    title('Real behavioral correlations');

    subplot(2,3,6); faFeat=nan(1,nf);
    for k=1:nf, notInj=~TM(:,k); faFeat(k)=mean(detMatch(notInj,k),'omitnan'); end
    bar(1:nf,faFeat,'FaceColor',[0.8 0.4 0.4]); yline(0.05,'k--','0.05');
    set(gca,'XTick',1:nf,'XTickLabel',feat,'XTickLabelRotation',25,'TickDir','out'); ylabel('False-alarm rate');
    ylim([0 max(0.2,max(faFeat)*1.2)]); title('Per-feature false positives'); box off;

    if useTtest, sigStr='t-test'; else, sigStr=sprintf('%d perms',nPerms); end
    if useNested, lamStr='\lambda=CV-tuned'; else, lamStr=sprintf('\\lambda=%.3g',lambdaVal); end
    sgtitle(sprintf('GLM ground-truth recovery | %s | %s | %s | green=injected',envChoice,lamStr,sigStr),'FontWeight','bold');
end

function cmap=redbluemap()
    n=64; half=n/2; r=[linspace(0,1,half),ones(1,half)]'; g=[linspace(0,1,half),linspace(1,0,half)]'; b=[ones(1,half),linspace(1,0,half)]'; cmap=[r g b];
end

function closeIfValid(d), if isvalid(d), close(d); end, end

function progressTick(d, nJobs, doReset)
    persistent cnt
    if isempty(cnt) || doReset, cnt = 0; end
    if doReset, return; end
    cnt = cnt + 1;
    if isvalid(d), d.Value = 0.15 + 0.8*(cnt/nJobs); d.Message = sprintf('Fitting %d/%d...', cnt, nJobs); end
end

function s = ternaryStr(cond,a,b)
    if cond, s=a; else, s=b; end
end
