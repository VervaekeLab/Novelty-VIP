function out = vipGLM_fitCell(sData, dff_all, cellID, envChoice, lambdaVal, nPerms, lambdaMode, sigTest, lagMin, lagMax, nLags, nPosRBF, posBasis)
% VIPGLM_FITCELL  Headless per-cell encoding-model (GLM) fit -- the core engine.
%
% AIM: fit the ridge-regression encoding model for ONE cell in ONE context and
% return its goodness-of-fit and per-feature unique variance. This is the
% computational heart of the toolkit; vipGLM_cohort and vipGLM_session call it
% on every cell. No figures, no GUI. Results match the interactive inspector.
%
% INPUT:  in-memory arguments (not a file prompt) -- a loaded sData struct and a
%         dff matrix; called by the orchestrators, not usually run by hand.
% OUTPUT: struct 'out' with full-model CV R^2, per-feature unique dR2,
%         significance, primary feature, and the held-out fit internals.
%
% Key options:
%   nPosRBF  : number of position basis functions (15 for RBF; 100 for Turi binary; 0 = drop position).
%   posBasis : 'binary' (Turi 2019 non-overlapping occupancy bins, default) | 'rbf' (Gaussian).
%
%   out = vipGLM_fitCell(sData, dff_all, cellID, envChoice, lambdaVal, nPerms, lambdaMode, sigTest, lagMin, lagMax, nLags)
%
%   Inputs
%     sData      - loaded session struct
%     dff_all    - nCells x nFrames dF/F matrix (sData.imdata.roiSignals(2).dff)
%     cellID     - row index into dff_all
%     envChoice  - 'Familiar Context' | 'Novel Context'
%     lambdaVal  - ridge penalty (used when lambdaMode='fixed')
%     nPerms     - circular-shift permutations for the null (sigTest='perm')
%     lambdaMode - 'fixed' (default) | 'nested'  (nested-CV grid search for lambda)
%     sigTest    - 'perm' (default) | 'ttest'    (paired t-test on per-fold MSE)
%
%   Output struct fields
%     R2_full          scalar pooled-CV R^2 of the full model
%     dR2_unique       1x6 unique variance (drop-one-out) per feature
%     single_R2        1x6 single-feature CV R^2
%     p_vals           1x6 significance p-values (perm or t-test)
%     feat_names       1x6 cellstr
%     sig              1x6 logical (p<0.05 & dR2>0.01)  -- the diversity mask
%     lambdaUsed       lambda actually used (= lambdaVal, or nested-CV choice)
%     nFramesUsed      frames in the analysis mask
%     Fs               sampling rate
%     hasPupil         logical, false if pupil absent (Pupil col then NaN'd)
%     ok               logical, false if cell/session could not be fit
%     msg              diagnostic string

    if nargin<7 || isempty(lambdaMode), lambdaMode='fixed'; end
    if nargin<8 || isempty(sigTest),    sigTest='perm';     end
    if nargin<9  || isempty(lagMin), lagMin=-3; end   % lag window start (s)
    if nargin<10 || isempty(lagMax), lagMax=3;  end   % lag window end (s)
    if nargin<11 || isempty(nLags),  nLags=5;   end   % number of lag points
    if nargin<12 || isempty(nPosRBF), nPosRBF=15; end  % position basis count (0 = no position)
    if nargin<13 || isempty(posBasis), posBasis='binary'; end  % 'binary' (Turi) | 'rbf'
    useNested = strcmpi(lambdaMode,'nested');
    useTtest  = strcmpi(sigTest,'ttest');

    out = struct('R2_full',NaN,'dR2_unique',nan(1,6),'single_R2',nan(1,6), ...
        'p_vals',nan(1,6),'feat_names',{{'Velocity','Acceleration','Licking', ...
        'Reward','Position','Pupil'}},'sig',false(1,6),'lambdaUsed',lambdaVal,'nFramesUsed',0, ...
        'Fs',NaN,'hasPupil',false,'pupilDetFrac',NaN,'ok',false,'msg','');

    try
        % ---- 1. DIMENSION BOUNDING ----
        dff_len = size(dff_all,2); time_len = length(sData.behavior.signals.timeDs);
        pos_len = length(sData.behavior.signals.corridorPositionDs);
        vel_len = length(sData.behavior.signals.velocityDs);
        N = min([dff_len, time_len, pos_len, vel_len]);

        time = sData.behavior.signals.timeDs(1:N); time = time(:);
        Fs = 1/mean(diff(time),'omitnan');
        if isnan(Fs) || Fs<=0, out.msg='bad Fs'; return; end
        out.Fs = Fs;

        dff  = double(dff_all(cellID,1:N)'); dff = dff(:);
        pos  = double(sData.behavior.signals.corridorPositionDs(1:N)); pos = pos(:);
        vel  = double(sData.behavior.signals.velocityDs(1:N)); vel = vel(:);
        lick = double(sData.behavior.signals.lickIfreqDs(1:N)); lick = lick(:);

        if length(vel)>2, vel(1:2)=vel(3); end
        acc = gradient(vel)*Fs; acc = acc(:);

        % ---- reward (downsampled vector preferred; verified DAQ fallback) ----
        valve = zeros(N,1);
        if isfield(sData.behavior,'signals') && isfield(sData.behavior.signals,'rewardDs')
            vt = double(sData.behavior.signals.rewardDs);
            valve(1:min(N,length(vt))) = vt(1:min(N,length(vt))) > 0;
        elseif isfield(sData,'daqdata') && isfield(sData.daqdata,'frameIndex') && isfield(sData.daqdata,'waterValve')
            vr = double(sData.daqdata.waterValve > 0.5);
            ft = find(diff(sData.daqdata.frameIndex)==1)+1;
            lim = min(N,length(ft)-1);
            for i=1:lim
                if any(vr(ft(i):ft(i+1)-1)), valve(i)=1; end
            end
            if N<=length(ft) && any(vr(ft(N):end)), valve(N)=1; end
        else
            out.msg='no reward data'; return;
        end
        valve = valve(:);

        % ---- pupil ----  (track which frames were ACTUALLY detected, pre-fill)
        pupil_detected = false(N,1);
        try
            if isfield(sData.behavior,'pupildata') && isfield(sData.behavior.pupildata,'diameter')
                pd = double(sData.behavior.pupildata.diameter); pd(isnan(pd))=inf;
                pfi = find(diff(sData.behavior.signals.pupilCamFrameIndexCorr)==1)+1;
                ifi = find(diff(sData.daqdata.frameIndex)==1);
                pUs = nan(size(sData.behavior.signals.pupilCamFrameIndexCorr));
                vl = min(length(pfi),length(pd));
                pUs(pfi(1:vl)) = pd(1:vl);
                detUs = isfinite(pUs);                 % detected before any interpolation
                pUs = fillmissing(pUs,'nearest'); pUs(isinf(pUs))=nan;
                pr = double(pUs(ifi)); dr = detUs(ifi);
                if length(pr)>=N, pupil=pr(1:N); det=dr(1:N);
                else, pupil=[pr(:);nan(N-length(pr),1)]; det=[dr(:);false(N-length(pr),1)]; end
                pupil = pupil(:); pupil_detected = det(:);
            else, pupil=nan(N,1); end
        catch, pupil=nan(N,1); end
        if all(isnan(pupil))
            pupil_glm=zeros(N,1); out.hasPupil=false;
        else
            pupil_glm=fillmissing(pupil,'nearest'); pupil_glm(isnan(pupil_glm))=mean(pupil_glm,'omitnan');
            out.hasPupil=true;   % provisional; re-checked within-context below
        end

        feat_names = out.feat_names;

        % ---- 3. CONTEXT MASK + CV CHUNKS ----
        ctx_raw = mode(sData.behavior.trialMatrices.contextMatrix,2);
        ctx = round(fillmissing(ctx_raw,'nearest')); ctx=ctx(:);
        trialIdx = sData.behavior.signals.trialIndexesDs(1:N); trialIdx=trialIdx(:);
        vtr = trialIdx>0 & trialIdx<=length(ctx);
        cctx = nan(N,1); cctx(vtr)=ctx(trialIdx(vtr)); cctx=fillmissing(cctx,'nearest');
        fam_id = ctx(1); uctx=unique(ctx); nov_id=uctx(uctx~=fam_id);

        if strcmp(envChoice,'Familiar Context')
            mask = (cctx==fam_id);   % all familiar frames (both blocks)
        else
            if isempty(nov_id), out.msg='no novel context'; return; end
            mask = (cctx==nov_id(1));
        end
        mask = mask(:);
        if sum(mask) < 50, out.msg='too few frames in context'; return; end
        out.nFramesUsed = sum(mask);

        % Record pupil detection fraction within the context as a DIAGNOSTIC only
        % (no longer used to gate hasPupil -- pupil data verified reliable). Anyone
        % can still filter on res.cellTable.pupilDetFrac post-hoc if desired.
        if out.hasPupil
            out.pupilDetFrac = mean(pupil_detected(mask));
        else
            out.pupilDetFrac = 0;
        end

        N_mask=sum(mask); k_folds=10; chunk_size=round(Fs*10);
        num_chunks=ceil(N_mask/chunk_size);
        fa=repmat(1:k_folds,1,ceil(num_chunks/k_folds)); cv_idx=zeros(N_mask,1);
        for i=1:num_chunks
            cv_idx(((i-1)*chunk_size+1):min(i*chunk_size,N_mask))=fa(i);
        end

        % ---- 4. DESIGN MATRIX (built ENTIRELY within the masked context) ----
        % Extract the active context first, then z-score / convolve / shift / RBF on
        % the masked data only. Nothing crosses the context boundary, so the real
        % regressors and the permutation null are constructed identically.
        zc = @(v) (v-mean(v,'omitnan'))./max(std(v,'omitnan'),1e-8);

        vel_m   = vel(mask);   acc_m = acc(mask);   lick_m = lick(mask);
        valve_m = valve(mask); pos_m = pos(mask);   pupil_m = pupil_glm(mask);
        dff_m   = dff(mask);
        Nm = sum(mask);

        dff_z = zc(dff_m); dff_z(isnan(dff_z))=0;

        sigma_frames=round(0.125*Fs); kernel=gausswin(round(sigma_frames*3)); kernel=kernel/sum(kernel);
        shifts=round(linspace(lagMin,lagMax,nLags)*Fs);   % lag window [lagMin,lagMax], nLags points

        X_vel  =make_shifted_matrix(conv(zc(vel_m),  kernel,'same'),shifts);
        X_acc  =make_shifted_matrix(conv(zc(acc_m),  kernel,'same'),shifts);
        X_lick =make_shifted_matrix(conv(zc(lick_m), kernel,'same'),shifts);
        X_rew  =make_shifted_matrix(conv(zc(valve_m),kernel,'same'),shifts);
        X_pupil=make_shifted_matrix(conv(zc(pupil_m),kernel,'same'),shifts);

        num_bases=nPosRBF;
        pos_edges=[];
        if num_bases>0
            if strcmpi(posBasis,'binary')
                % Turi 2019: divide the track into num_bases non-overlapping bins,
                % each a binary indicator (1 when the animal occupies that bin).
                pos_edges=linspace(min(pos),max(pos),num_bases+1);
                pos_centers=(pos_edges(1:end-1)+pos_edges(2:end))/2; pos_width=NaN;
                pb=discretize(pos_m,pos_edges);   % NaN (out-of-range/dropped) -> all-zero pos row below
                X_pos=zeros(Nm,num_bases);
                for b=1:num_bases, X_pos(:,b)=double(pb==b); end
            else
                % Gaussian RBF basis (smooth, overlapping) -- alternative to Turi.
                pos_centers=linspace(0,max(pos),num_bases);
                pos_width=(max(pos)/num_bases)*1.2;
                X_pos=zeros(Nm,num_bases);
                for b=1:num_bases, X_pos(:,b)=exp(-((pos_m-pos_centers(b)).^2)/(2*pos_width^2)); end
            end
        else
            pos_centers=[]; pos_width=1; X_pos=zeros(Nm,0);   % position removed from model
        end

        % Standardize the position basis so the SINGLE ridge penalty falls on it
        % the same way it falls on the (already z-scored) continuous features.
        % Binary occupancy columns sit at SD~0.1; without this, lambda crushes
        % position and biases its unique dR2 downward. Applied IDENTICALLY to the
        % permutation-null position block below: circular shift preserves each
        % column's mean/SD, so null and real scaling coincide and the test stays
        % valid. All-zero columns (unvisited bins) -> the 1e-8 guard leaves them 0.
        if size(X_pos,2)>0
            X_pos = (X_pos - mean(X_pos,1))./max(std(X_pos,0,1),1e-8);
        end

        % Robust block indexing via cumulative widths -- safe when any feature
        % block is 0-width (e.g. nPosRBF=0 removes position). The old max(idx)+1
        % chain breaks here: max([]) returns [] and corrupts every later block.
        cols  = [size(X_vel,2), size(X_acc,2), size(X_lick,2), size(X_rew,2), size(X_pos,2), size(X_pupil,2)];
        ends  = cumsum(cols);
        starts= [1, ends(1:end-1)+1];
        idx_vel  = starts(1):ends(1);
        idx_acc  = starts(2):ends(2);
        idx_lick = starts(3):ends(3);
        idx_rew  = starts(4):ends(4);
        idx_pos  = starts(5):ends(5);
        idx_pupil= starts(6):ends(6);

        Xm=[X_vel,X_acc,X_lick,X_rew,X_pos,X_pupil]; Xm(isnan(Xm))=0;
        feat_indices={idx_vel,idx_acc,idx_lick,idx_rew,idx_pos,idx_pupil};

        % masked raw behaviors for the within-context permutation null
        raw_masked = {vel_m, acc_m, lick_m, valve_m, pos_m, pupil_m};

        % ---- 5. FIT + PARTITION + SIGNIFICANCE ----
        ym=dff_z;

        % resolve lambda: fixed value, or CV-tuned grid search on this cell
        % (NB: tuned on the same folds used for evaluation -- not true nested CV)
        if useNested
            grid=[0.01 0.1 0.5 1 2 5 10 50]; gr=zeros(numel(grid),1);
            for gi=1:numel(grid), gr(gi)=fit_cv_ridge_pooled(Xm,ym,cv_idx,k_folds,grid(gi)); end
            [~,bi]=max(gr); lambdaUsed=grid(bi);
        else
            lambdaUsed=lambdaVal;
        end
        out.lambdaUsed=lambdaUsed;

        R2_full = fit_cv_ridge_pooled(Xm,ym,cv_idx,k_folds,lambdaUsed);
        out.R2_full = R2_full;

        nf=6; dR2=zeros(1,nf); sR2=zeros(1,nf); pv=ones(1,nf);
        if R2_full<=0
            out.dR2_unique=dR2; out.single_R2=sR2; out.p_vals=pv;
            out.sig=false(1,nf); out.ok=true; out.msg='unresponsive'; 
            if ~out.hasPupil, out.dR2_unique(6)=NaN; out.single_R2(6)=NaN; out.p_vals(6)=NaN; end
            return;
        end

        for fI=1:nf
            if isempty(feat_indices{fI})
                dR2(fI)=0; sR2(fI)=0; pv(fI)=1;   % feature removed from model
                continue;
            end
            keep=setdiff(1:size(Xm,2),feat_indices{fI});
            R2_red=fit_cv_ridge_pooled(Xm(:,keep),ym,cv_idx,k_folds,lambdaUsed);
            dR2(fI)=R2_full-R2_red;

            if useTtest
                % paired t-test on per-fold MSE (full vs reduced), no permutations
                [mse_red,mse_full]=fold_mse_compare(Xm(:,keep),Xm,ym,cv_idx,k_folds,lambdaUsed);
                try [~,pv(fI)]=ttest(mse_red,mse_full,'Tail','right'); catch, pv(fI)=1; end
            else
                ndr=zeros(nPerms,1); rv=raw_masked{fI};   % within-context raw vector
                isPos=strcmp(feat_names{fI},'Position');
                for p=1:nPerms
                    sa=randi([round(Nm*0.1),round(Nm*0.9)]);
                    rs=circshift(rv,sa,1);                 % circular shift WITHIN the context
                    if ~isPos
                        Xb=make_shifted_matrix(conv(zc(rs),kernel,'same'),shifts);
                    elseif strcmpi(posBasis,'binary')
                        pbb=discretize(rs,pos_edges);   % NaN -> all-zero pos row (matches real matrix)
                        Xb=zeros(Nm,num_bases);
                        for b=1:num_bases, Xb(:,b)=double(pbb==b); end
                    else
                        Xb=zeros(Nm,num_bases);
                        for b=1:num_bases, Xb(:,b)=exp(-((rs-pos_centers(b)).^2)/(2*pos_width^2)); end
                    end
                    % Mirror the real-matrix position standardization (see above).
                    % Continuous features are already z-scored via zc(rs); position
                    % must be too, or the null is built differently from the real X.
                    if isPos && size(Xb,2)>0
                        Xb = (Xb - mean(Xb,1))./max(std(Xb,0,1),1e-8);
                    end
                    Xs=Xm; Xs(:,feat_indices{fI})=Xb; Xs(isnan(Xs))=0;
                    ndr(p)=fit_cv_ridge_pooled(Xs,ym,cv_idx,k_folds,lambdaUsed)-R2_red;
                end
                pv(fI)=(sum(ndr>=dR2(fI))+1)/(nPerms+1);
            end

            sR2(fI)=fit_cv_ridge_pooled(Xm(:,feat_indices{fI}),ym,cv_idx,k_folds,lambdaUsed);
        end

        out.dR2_unique=dR2; out.single_R2=sR2; out.p_vals=pv;
        out.sig = (pv<0.05) & (dR2>0.01);

        % If pupil absent, mark its slots NaN so it isn't counted as a real predictor
        if ~out.hasPupil
            out.dR2_unique(6)=NaN; out.single_R2(6)=NaN; out.p_vals(6)=NaN; out.sig(6)=false;
        end
        out.ok=true; out.msg='ok';

    catch ME
        out.ok=false; out.msg=ME.message;
    end
end

% ===== helpers (identical math to the inspector) =====
function X_shifted = make_shifted_matrix(vec, shifts)
    N=length(vec); X_shifted=zeros(N,length(shifts));
    for i=1:length(shifts)
        s=shifts(i);
        if s==0, X_shifted(:,i)=vec;
        elseif s>0, X_shifted(:,i)=[zeros(s,1);vec(1:end-s)];
        else, sa=abs(s); X_shifted(:,i)=[vec(sa+1:end);zeros(sa,1)]; end
    end
end

function R2 = fit_cv_ridge_pooled(X, y, cv_idx, k_folds, lambdaVal)
    N=length(y); yp=nan(N,1);
    for f=1:k_folds
        te=(cv_idx==f); tr=(cv_idx~=f & cv_idx>0);
        if sum(tr)==0||sum(te)==0, continue; end
        Xtr=[ones(sum(tr),1),X(tr,:)]; Xte=[ones(sum(te),1),X(te,:)];
        pen=lambdaVal*eye(size(Xtr,2)); pen(1,1)=0;
        b=(Xtr'*Xtr+pen)\(Xtr'*y(tr));
        yp(te)=Xte*b;
    end
    vi=~isnan(yp)&cv_idx>0;
    if sum(vi)==0, R2=0; return; end
    sst=sum((y(vi)-mean(y(vi))).^2); sse=sum((y(vi)-yp(vi)).^2);
    if sst>0, R2=1-(sse/sst); else, R2=0; end
end

function [mse_red, mse_full] = fold_mse_compare(Xred, Xfull, y, cv_idx, k_folds, lambdaVal)
% Per-fold test MSE for reduced and full models (for the t-test significance mode)
    mse_red=nan(k_folds,1); mse_full=nan(k_folds,1);
    for f=1:k_folds
        te=(cv_idx==f); tr=(cv_idx~=f & cv_idx>0);
        if sum(tr)==0||sum(te)==0, continue; end
        Xtr_r=[ones(sum(tr),1),Xred(tr,:)]; Xte_r=[ones(sum(te),1),Xred(te,:)];
        pr=lambdaVal*eye(size(Xtr_r,2)); pr(1,1)=0; br=(Xtr_r'*Xtr_r+pr)\(Xtr_r'*y(tr));
        Xtr_f=[ones(sum(tr),1),Xfull(tr,:)]; Xte_f=[ones(sum(te),1),Xfull(te,:)];
        pf=lambdaVal*eye(size(Xtr_f,2)); pf(1,1)=0; bf=(Xtr_f'*Xtr_f+pf)\(Xtr_f'*y(tr));
        mse_red(f)=mean((y(te)-Xte_r*br).^2);
        mse_full(f)=mean((y(te)-Xte_f*bf).^2);
    end
    mse_red=mse_red(~isnan(mse_red)); mse_full=mse_full(~isnan(mse_full));
end
