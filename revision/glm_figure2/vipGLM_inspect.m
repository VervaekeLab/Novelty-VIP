function vipGLM_inspect()
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_INSPECT  Interactive single-cell encoding-model inspector (GUI).
%
% AIM: explore one cell at a time -- pick a session and cell, see its GLM fit,
% per-feature unique dR2, and tuning. The original interactive tool the headless
% engine (vipGLM_fitCell) was derived from; useful for eyeballing individual cells.
%
% INPUT:  ONE SESSION DATA file (sData .mat), via prompt.
% OUTPUT: interactive inspector figure (no files saved).
    % 1. File Selection
    userHome = getenv('HOME'); 
    folderPath = fullfile(userHome, 'Desktop', 'VIP-GCaMP-sData');
    if ~exist(folderPath, 'dir'), folderPath = fullfile(userHome, 'Desktop'); end
    
        uiwait(msgbox({'Select ONE SESSION DATA file.','','A single m8xxx session .mat file.'},'Session data (one file)','help'));
    [fileName, pathName] = uigetfile('*.mat','ONE SESSION DATA file (m8xxx .mat)',dfSess);
    if isequal(fileName, 0), return; end
    
    fprintf('Loading %s...\n', fileName);
    data = load(fullfile(pathName, fileName));
    sData = data.sData;
    
    if isfield(sData.imdata, 'roiSignals') && length(sData.imdata.roiSignals) >= 2
        dff_all = sData.imdata.roiSignals(2).dff;
    else
        error('Could not find VIP neural data in sData.imdata.roiSignals(2).dff');
    end
    nCells = size(dff_all, 1);
    
    % 2. Setup GUI 
    uiFig = uifigure('Name', ['GLM + Spatial Inspector: ' fileName], 'Position', [100 100 320 510]);
    envDropdown = uidropdown(uiFig, 'Position', [20 450 280 30], 'Items', {'Familiar Context', 'Novel Context'});
    spinner = uispinner(uiFig, 'Position', [20 405 280 30], 'Limits', [1 nCells], 'Value', 1);
    statsDropdown = uidropdown(uiFig, 'Position', [20 360 280 30], 'Items', {'Exploratory Mode (Fast t-test)', 'Publication Mode (500 Permutations)'});
    lambdaModeDropdown = uidropdown(uiFig, 'Position', [20 315 280 30], 'Items', {'Fixed lambda', 'Nested CV lambda'});
    parDropdown = uidropdown(uiFig, 'Position', [20 270 280 30], 'Items', {'Parallel: Yes', 'Parallel: No'}, 'Value', 'Parallel: Yes');
    uilabel(uiFig, 'Position', [20 225 130 30], 'Text', 'Ridge \lambda:');
    lambdaField = uieditfield(uiFig, 'numeric', 'Position', [150 225 150 30], 'Value', 0.5);
    uilabel(uiFig, 'Position', [20 185 130 30], 'Text', 'Lag min (s):');
    lagMinField = uieditfield(uiFig, 'numeric', 'Position', [150 185 150 30], 'Value', -3, 'Limits', [-10 10]);
    uilabel(uiFig, 'Position', [20 145 130 30], 'Text', 'Lag max (s):');
    lagMaxField = uieditfield(uiFig, 'numeric', 'Position', [150 145 150 30], 'Value', 3, 'Limits', [-10 10]);
    uilabel(uiFig, 'Position', [20 105 130 30], 'Text', '# lag points:');
    nLagField = uieditfield(uiFig, 'numeric', 'Position', [150 105 150 30], 'Value', 5, 'Limits', [1 21], 'RoundFractionalValues', 'on');
    uibutton(uiFig, 'Text', 'Fit Model & Plot Selected Cell', 'Position', [20 20 280 55], ...
        'ButtonPushedFcn', @(src, event) plotDualInspector(sData, dff_all, spinner.Value, envDropdown.Value, ...
        statsDropdown.Value, lambdaField.Value, lambdaModeDropdown.Value, parDropdown.Value, ...
        lagMinField.Value, lagMaxField.Value, round(nLagField.Value), uiFig));
end

function plotDualInspector(sData, dff_all, cellID, envChoice, statsMode, lambdaVal, lambdaMode, parMode, lagMin, lagMax, nLags, uiFig)
    d = uiprogressdlg(uiFig, 'Title', sprintf('Analyzing Cell #%d', cellID), 'Message', 'Extracting and aligning signals...', 'Value', 0.05);

    useNestedCV = contains(lambdaMode, 'Nested');
    useParallel = contains(parMode, 'Yes');

    try
        % 1. ROBUST DIMENSION BOUNDING
        dff_len = size(dff_all, 2); time_len = length(sData.behavior.signals.timeDs);
        pos_len = length(sData.behavior.signals.corridorPositionDs); vel_len = length(sData.behavior.signals.velocityDs);
        N = min([dff_len, time_len, pos_len, vel_len]);
        
        time = sData.behavior.signals.timeDs(1:N); time = time(:);
        Fs = 1 / mean(diff(time), 'omitnan'); 
        if isnan(Fs) || Fs <= 0, error('Could not determine valid sampling rate (Fs) from time array.'); end
        
        dff = double(dff_all(cellID, 1:N)'); dff = dff(:);
        pos = double(sData.behavior.signals.corridorPositionDs(1:N)); pos = pos(:);
        vel = double(sData.behavior.signals.velocityDs(1:N)); vel = vel(:);
        lick = double(sData.behavior.signals.lickIfreqDs(1:N)); lick = lick(:);
        
        % Patch initial derivative artifact and scale to proper physical units (cm/s^2)
        if length(vel) > 2, vel(1:2) = vel(3); end
        acc = gradient(vel) * Fs; acc = acc(:);
        
        % Strict requirement for downsampled reward vector to prevent misalignment
        valve = zeros(N, 1);
        if isfield(sData.behavior, 'signals') && isfield(sData.behavior.signals, 'rewardDs')
            valve_tmp = double(sData.behavior.signals.rewardDs);
            valve(1:min(N, length(valve_tmp))) = valve_tmp(1:min(N, length(valve_tmp))) > 0;
        elseif isfield(sData, 'daqdata') && isfield(sData.daqdata, 'frameIndex') && isfield(sData.daqdata, 'waterValve')
            fprintf('Notice: rewardDs missing. Executing verified sub-frame DAQ integration...\n');
            valve_raw = double(sData.daqdata.waterValve > 0.5); 
            f_trans = find(diff(sData.daqdata.frameIndex) == 1) + 1; 
            limit = min(N, length(f_trans)-1);
            for i = 1:limit
                if any(valve_raw(f_trans(i) : f_trans(i+1)-1)), valve(i) = 1; end
            end
            if N <= length(f_trans) && any(valve_raw(f_trans(N) : end)), valve(N) = 1; end
        else
            error('No reward data found in behavior.signals or daqdata.');
        end
        valve = valve(:);
        
        try
            if isfield(sData.behavior, 'pupildata') && isfield(sData.behavior.pupildata, 'diameter')
                pupilDiameter = double(sData.behavior.pupildata.diameter); pupilDiameter(isnan(pupilDiameter)) = inf; 
                pupilFrameInd = find(diff(sData.behavior.signals.pupilCamFrameIndexCorr)==1)+1; 
                imagingFrameInd = find(diff(sData.daqdata.frameIndex)==1);
                pupilDiameterUs = nan(size(sData.behavior.signals.pupilCamFrameIndexCorr));
                valid_len = min(length(pupilFrameInd), length(pupilDiameter));
                pupilDiameterUs(pupilFrameInd(1:valid_len)) = pupilDiameter(1:valid_len);
                pupilDiameterUs = fillmissing(pupilDiameterUs, 'nearest'); pupilDiameterUs(isinf(pupilDiameterUs)) = nan;
                pupil_raw = double(pupilDiameterUs(imagingFrameInd));
                if length(pupil_raw) >= N, pupil = pupil_raw(1:N); else, pupil = [pupil_raw(:); nan(N - length(pupil_raw), 1)]; end
                pupil = pupil(:);
            else, pupil = nan(N, 1); end
        catch, pupil = nan(N, 1); end
        if all(isnan(pupil)), pupil_glm = zeros(N, 1); else, pupil_glm = fillmissing(pupil, 'nearest'); pupil_glm(isnan(pupil_glm)) = mean(pupil_glm, 'omitnan'); end
        
        feat_names = {'Velocity', 'Acceleration', 'Licking', 'Reward', 'Position', 'Pupil'};
        
        % 3. CONTEXT MASKING & TIME CHUNKING
        ctx_raw = mode(sData.behavior.trialMatrices.contextMatrix, 2); ctx = round(fillmissing(ctx_raw, 'nearest')); ctx = ctx(:);
        trialIdx = sData.behavior.signals.trialIndexesDs(1:N); trialIdx = trialIdx(:);
        valid_trials = trialIdx > 0 & trialIdx <= length(ctx);
        continuous_ctx = nan(N, 1); continuous_ctx(valid_trials) = ctx(trialIdx(valid_trials)); continuous_ctx = fillmissing(continuous_ctx, 'nearest');
        fam_id = ctx(1); unique_ctx = unique(ctx); nov_id = unique_ctx(unique_ctx ~= fam_id);
        
        if strcmp(envChoice, 'Familiar Context')
            mask = (continuous_ctx == fam_id);   % all familiar frames (both blocks)
        else
            if isempty(nov_id), error('No Novel Context found in this session.'); end
            mask = (continuous_ctx == nov_id(1));
        end
        mask = mask(:); 
        
        N_mask = sum(mask); k_folds = 10; chunk_size = round(Fs * 10); num_chunks = ceil(N_mask / chunk_size);
        fold_assignments = repmat(1:k_folds, 1, ceil(num_chunks/k_folds)); cv_idx = zeros(N_mask, 1);
        for i = 1:num_chunks, cv_idx(((i-1)*chunk_size + 1) : min(i*chunk_size, N_mask)) = fold_assignments(i); end
        
        d.Value = 0.10; d.Message = 'Building Design Matrix...';
    
        % 4. BUILD GLM DESIGN MATRIX  (built ENTIRELY within the masked context)
        % Mask raw behavior first, then z-score / convolve / shift / RBF on masked
        % data only, so nothing crosses the context boundary and the permutation
        % null is constructed identically to the real regressors.
        ctx_z = @(v, m) (v - mean(v(m), 'omitnan')) ./ max(std(v(m), 'omitnan'), 1e-8); % (kept for plots)
        zc = @(v) (v - mean(v,'omitnan')) ./ max(std(v,'omitnan'),1e-8);
        dff_z = ctx_z(dff, mask);   % full-length z (for the trace plot)

        vel_m=vel(mask); acc_m=acc(mask); lick_m=lick(mask);
        valve_m=valve(mask); pos_m=pos(mask); pupil_m=pupil_glm(mask);
        Nm=sum(mask);
        ym_fit = zc(dff(mask)); ym_fit(isnan(ym_fit))=0;
        raw_masked = {vel_m, acc_m, lick_m, valve_m, pos_m, pupil_m};

        sigma_frames = round(0.125 * Fs); kernel = gausswin(round(sigma_frames*3)); kernel = kernel/sum(kernel);
        
        % +/- 3 s symmetric window (matches Turi et al. 2019)
        shifts = round(linspace(lagMin, lagMax, nLags) * Fs);   % lag window [lagMin,lagMax], nLags points
        
        X_vel   = make_shifted_matrix(conv(zc(vel_m),   kernel, 'same'), shifts);
        X_acc   = make_shifted_matrix(conv(zc(acc_m),   kernel, 'same'), shifts);
        X_lick  = make_shifted_matrix(conv(zc(lick_m),  kernel, 'same'), shifts);
        X_rew   = make_shifted_matrix(conv(zc(valve_m), kernel, 'same'), shifts);
        X_pupil = make_shifted_matrix(conv(zc(pupil_m), kernel, 'same'), shifts);
        
        num_bases = 15; pos_centers = linspace(0, max(pos), num_bases); pos_width = (max(pos) / num_bases) * 1.2; 
        X_pos = zeros(Nm, num_bases);
        for b = 1:num_bases, X_pos(:, b) = exp(-((pos_m - pos_centers(b)).^2) / (2 * pos_width^2)); end
        
        idx_vel   = 1:size(X_vel,2);
        idx_acc   = max(idx_vel)+1 : max(idx_vel)+size(X_acc,2);
        idx_lick  = max(idx_acc)+1 : max(idx_acc)+size(X_lick,2);
        idx_rew   = max(idx_lick)+1 : max(idx_lick)+size(X_rew,2);
        idx_pos   = max(idx_rew)+1 : max(idx_rew)+num_bases;
        idx_pupil = max(idx_pos)+1 : max(idx_pos)+size(X_pupil,2);
        
        Xm_fit = [X_vel, X_acc, X_lick, X_rew, X_pos, X_pupil];
        Xm_fit(isnan(Xm_fit)) = 0;
        
        feat_indices = {idx_vel, idx_acc, idx_lick, idx_rew, idx_pos, idx_pupil};
        
        % Resolve lambda: fixed, or CV-tuned grid search (tuned on same folds -- not true nested)
        if useNestedCV
            d.Value = 0.13; d.Message = 'Selecting \lambda by CV...';
            lambdaUsed = select_lambda_nested(Xm_fit, ym_fit, cv_idx, k_folds);
        else
            lambdaUsed = lambdaVal;
        end
        
        d.Value = 0.15; d.Message = sprintf('Fitting GLM (\\lambda=%.3g)...', lambdaUsed);
    
        % 5. RUN GLM  (all inputs are within-context / masked length)
        y_target = dff_z; y_target(isnan(y_target)) = 0;   % full-length, for plotting
        [R2_full, y_pred_trace, dR2_unique, single_R2, p_vals, selectivity_label] = run_statistical_glm(...
            Xm_fit, ym_fit, feat_indices, feat_names, d, statsMode, cv_idx, k_folds, lambdaUsed, ...
            raw_masked, zc, kernel, shifts, num_bases, pos_centers, pos_width, Nm, useParallel);
        
        residual_variance = R2_full - sum(dR2_unique);
        
        vel_corr = corr(vel(mask), y_target(mask), 'rows', 'complete'); if vel_corr < 0, vel_tag = '(-)'; else, vel_tag = '(+)'; end
        if all(isnan(pupil_glm)), pupil_tag = ''; else
            pupil_corr = corr(pupil_glm(mask), y_target(mask), 'rows', 'complete');
            if pupil_corr < 0, pupil_tag = '(-)'; else, pupil_tag = '(+)'; end
        end
        feat_names_display = feat_names; feat_names_display{1} = sprintf('Velocity %s', vel_tag); feat_names_display{6} = sprintf('Pupil %s', pupil_tag);
        
        d.Value = 0.90; d.Message = 'Plotting...';
    
        % 6. CALCULATE SPATIAL TRIAL AVERAGES
        numBins = 50; binEdges = linspace(0, max(pos), numBins + 1); binCenters = linspace(0, max(pos), numBins);
        mDFF=nan(numBins,1); sDFF=nan(numBins,1); mVel=nan(numBins,1); sVel=nan(numBins,1); mAcc=nan(numBins,1); sAcc=nan(numBins,1); mLick=nan(numBins,1); sLick=nan(numBins,1); mValve=nan(numBins,1); mPupil=nan(numBins,1); sPupil=nan(numBins,1);
        for b = 1:numBins
            idx = (pos >= binEdges(b) & pos < binEdges(b+1)) & mask; n_pts = sum(idx);
            if n_pts > 0
                mDFF(b) = mean(dff_z(idx), 'omitnan'); sDFF(b) = std(dff_z(idx), 'omitnan')/sqrt(n_pts); 
                mVel(b) = mean(vel(idx), 'omitnan'); sVel(b) = std(vel(idx), 'omitnan')/sqrt(n_pts);
                mAcc(b) = mean(acc(idx), 'omitnan'); sAcc(b) = std(acc(idx), 'omitnan')/sqrt(n_pts); 
                mLick(b) = mean(lick(idx), 'omitnan'); sLick(b) = std(lick(idx), 'omitnan')/sqrt(n_pts);
                mValve(b) = mean(valve(idx), 'omitnan'); mPupil(b) = mean(pupil(idx), 'omitnan'); sPupil(b) = std(pupil(idx), 'omitnan')/sqrt(n_pts);
            end
        end
        
        % 7. PLOTTING
        f = figure(888); clf(f); set(f, 'Color', 'w', 'Name', sprintf('GLM Diagnostic: Cell #%d', cellID), 'Position', [50 50 1400 850]);
        
        t_mask = (1:sum(mask))' / Fs; ax_trace = subplot(4, 2, [1 3 5]); hold on;
        plot(t_mask, y_target(mask), 'k', 'LineWidth', 1); plot(t_mask, y_pred_trace, 'r', 'LineWidth', 1.2);
        ylabel('Z-Scored \DeltaF/F'); title(sprintf('GLM Fit: Cell #%d (%s, \\lambda=%.3g)\n%s (Pooled Block R^2=%.2f)', cellID, envChoice, lambdaUsed, selectivity_label, R2_full));
        legend('Actual \DeltaF/F', 'Full GLM Prediction', 'Location', 'northeast'); axis tight; set(gca, 'TickDir', 'out');
        
        window_size = 60; set(ax_trace, 'XLim', [0, window_size]); slider_max = max(1, t_mask(end) - window_size);
        uicontrol('Style', 'slider', 'Min', 0, 'Max', slider_max, 'Value', 0, 'Units', 'normalized', 'Position', [0.13 0.35 0.33 0.02], 'Callback', @(s,e) set(ax_trace, 'XLim', [s.Value, s.Value + window_size]));
        
        subplot(4, 2, 7); hold on;
        b = bar(1:6, [dR2_unique', single_R2'], 'grouped', 'EdgeColor', 'k');
        b(1).FaceColor = [0.2 0.6 0.8]; b(2).FaceColor = [0.7 0.8 0.9]; 
        for i = 1:6
            if p_vals(i) < 0.05 && dR2_unique(i) > 0.01
                x_pos = b(1).XEndPoints(i); text(x_pos, dR2_unique(i) + max([dR2_unique, single_R2])*0.05, '*', 'FontSize', 20, 'HorizontalAlignment', 'center', 'Color', 'r');
            end
        end
        set(gca, 'XTick', 1:6, 'XTickLabel', feat_names_display, 'TickDir', 'out', 'XTickLabelRotation', 15);
        ylabel('Explained Variance (R^2)'); title(sprintf('Variance Partitioning (Residual / Shared: %.2f)', residual_variance)); box off;
        legend({'Unique (Drop-One-Out)', 'Total (Single-Feature Only)'}, 'Location', 'northwest');
        
        rewardLines = [100, 220]; 
        ax1 = subplot(4, 2, 2); hold on; plotShaded(binCenters, mDFF, sDFF, [0.3 0.3 0.3]); plot(binCenters, mDFF, 'k', 'LineWidth', 2); xline(rewardLines, '--r', 'Alpha', 0.5); ylabel('Z-Scored \DeltaF/F'); title('Continuous \DeltaF/F Average'); axis tight; set(gca, 'XTickLabel', []);
        ax2 = subplot(4, 2, 4); hold on; yyaxis left; plotShaded(binCenters, mVel, sVel, [0.2 0.5 0.9]); plot(binCenters, mVel, '-b', 'LineWidth', 2); ylabel('Vel (cm/s)'); yyaxis right; plotShaded(binCenters, mAcc, sAcc, [0.9 0.4 0.1]); plot(binCenters, mAcc, '-', 'Color', [0.8 0.3 0], 'LineWidth', 1.5); ylabel('Accel'); xline(rewardLines, '--k', 'Alpha', 0.3); title('Locomotion'); axis tight; set(gca, 'XTickLabel', []);
        ax3 = subplot(4, 2, 6); hold on; yyaxis left; plotShaded(binCenters, mLick, sLick, [0.8 0.2 0.8]); plot(binCenters, mLick, '-m', 'LineWidth', 2); ylabel('Lick Rate'); yyaxis right; bar(binCenters, mValve, 'FaceColor', [0.2 0.8 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.2); ylabel('Reward Dens.'); xline(rewardLines, '--k', 'Alpha', 0.3); title('Consummatory Behavior'); axis tight; set(gca, 'XTickLabel', []);
        ax4 = subplot(4, 2, 8); hold on; if all(isnan(mPupil)), text(mean(binCenters), 0, 'No Pupil', 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]); set(gca, 'YTick', []); else, plotShaded(binCenters, mPupil, sPupil, [0.2 0.7 0.7]); plot(binCenters, mPupil, '-', 'Color', [0 0.6 0.6], 'LineWidth', 2); ylabel('Pupil'); end; xlabel('Position (cm)'); xline(rewardLines, '--k', 'Alpha', 0.3); title('Arousal'); axis tight;
        linkaxes([ax1, ax2, ax3, ax4], 'x'); 
        
        if isvalid(d), close(d); end

    catch ME
        if isvalid(d), close(d); end
        uialert(uiFig, sprintf('Error: %s\n\nLine: %d', ME.message, ME.stack(1).line), 'Execution Error');
        fprintf('Error in plotDualInspector:\n%s\n', getReport(ME));
    end
end

% --- HELPERS ---
function plotShaded(x, y, e, col)
    v = ~isnan(y); 
    if any(v)
        fill([x(v), fliplr(x(v))], [y(v)'+e(v)', fliplr(y(v)'-e(v)')], col, 'EdgeColor', 'none', 'FaceAlpha', 0.3); 
    end
end

function X_shifted = make_shifted_matrix(vec, shifts)
    N = length(vec); X_shifted = zeros(N, length(shifts));
    for i = 1:length(shifts)
        s = shifts(i); if s == 0, X_shifted(:, i) = vec; elseif s > 0, X_shifted(:, i) = [zeros(s, 1); vec(1:end-s)]; else, s_abs = abs(s); X_shifted(:, i) = [vec(s_abs+1:end); zeros(s_abs, 1)]; end
    end
end

% Select lambda by nested CV: inner search over a grid, pick best mean test R^2.
function bestLambda = select_lambda_nested(X, y, cv_idx, k_folds)
    grid = [0.01 0.1 0.5 1 2 5 10 50];
    r2 = zeros(numel(grid),1);
    for gi = 1:numel(grid)
        r2(gi) = fit_cv_ridge_pooled(X, y, cv_idx, k_folds, grid(gi));
    end
    [~, bi] = max(r2);
    bestLambda = grid(bi);
    fprintf('Nested CV selected lambda = %.3g (pooled R^2=%.3f)\n', bestLambda, r2(bi));
end

function [R2_full, y_pred_full, delta_R2_unique, single_R2, p_vals, selectivity_label] = run_statistical_glm(X, y, feat_indices, feat_names, d, statsMode, cv_idx, k_folds, lambdaVal, raw_masked, zc, kernel, shifts, num_bases, pos_centers, pos_width, Nm, useParallel)
    isRigorous = contains(statsMode, 'Publication');
    
    [R2_full, y_pred_full] = fit_cv_ridge_pooled(X, y, cv_idx, k_folds, lambdaVal);
    
    n_feats = length(feat_names); delta_R2_unique = zeros(1, n_feats); single_R2 = zeros(1, n_feats); p_vals = ones(1, n_feats);
    if R2_full <= 0, selectivity_label = 'Unresponsive (R^2 \leq 0)'; return; end
    
    for f = 1:n_feats
        if isvalid(d), d.Value = 0.15 + (0.7 * (f / n_feats)); d.Message = sprintf('Evaluating feature %d of %d: %s...', f, n_feats, feat_names{f}); end
        cols_to_keep = setdiff(1:size(X,2), feat_indices{f});
        
        if isRigorous
            [R2_red, ~] = fit_cv_ridge_pooled(X(:, cols_to_keep), y, cv_idx, k_folds, lambdaVal);
            delta_R2_unique(f) = R2_full - R2_red; 
            
            n_perms = 500; null_dR2 = zeros(n_perms, 1);
            raw_vec = raw_masked{f};        % within-context raw vector
            feat_idx = feat_indices{f};
            is_pos = strcmp(feat_names{f}, 'Position');
            
            if useParallel
                parfor p = 1:n_perms
                    null_dR2(p) = one_perm(X, y, raw_vec, feat_idx, is_pos, zc, kernel, shifts, ...
                        num_bases, pos_centers, pos_width, Nm, cv_idx, k_folds, lambdaVal, R2_red);
                end
            else
                for p = 1:n_perms
                    null_dR2(p) = one_perm(X, y, raw_vec, feat_idx, is_pos, zc, kernel, shifts, ...
                        num_bases, pos_centers, pos_width, Nm, cv_idx, k_folds, lambdaVal, R2_red);
                end
            end
            p_vals(f) = (sum(null_dR2 >= delta_R2_unique(f)) + 1) / (n_perms + 1);
        else
            [R2_red, ~, fold_mse_red, fold_mse_full] = fit_cv_ridge_compare(X(:, cols_to_keep), X, y, cv_idx, k_folds, lambdaVal);
            delta_R2_unique(f) = R2_full - R2_red; 
            try [~, p_vals(f)] = ttest(fold_mse_red, fold_mse_full, 'Tail', 'right'); catch, p_vals(f) = 1; end
        end
        
        [R2_sing, ~] = fit_cv_ridge_pooled(X(:, feat_indices{f}), y, cv_idx, k_folds, lambdaVal);
        single_R2(f) = R2_sing;
    end
    
    valid_idx = find(p_vals < 0.05 & delta_R2_unique > 0.01);
    if isempty(valid_idx)
        if R2_full > 0.05, selectivity_label = sprintf('Correlated Coding (R^2=%.2f, no unique)', R2_full); else, selectivity_label = 'No Significant Unique Tuning'; end
    elseif length(valid_idx) == 1, selectivity_label = sprintf('Unique Contributor: %s', feat_names{valid_idx});
    else, selectivity_label = sprintf('Unique Contributors: %s', strjoin(feat_names(valid_idx), ', ')); end
end

% One within-context circular-shift permutation (shared by parallel and serial paths).
function ndr = one_perm(X, y, raw_vec, feat_idx, is_pos, zc, kernel, shifts, num_bases, pos_centers, pos_width, Nm, cv_idx, k_folds, lambdaVal, R2_red)
    shift_amt = randi([round(Nm*0.1), round(Nm*0.9)]); 
    raw_shuff = circshift(raw_vec, shift_amt, 1);     % shift WITHIN the context
    if ~is_pos
        X_block_shuff = make_shifted_matrix(conv(zc(raw_shuff), kernel, 'same'), shifts);
    else 
        X_block_shuff = zeros(Nm, num_bases);
        for b = 1:num_bases
            X_block_shuff(:, b) = exp(-((raw_shuff - pos_centers(b)).^2) / (2 * pos_width^2)); 
        end
    end
    X_shuff = X; X_shuff(:, feat_idx) = X_block_shuff; X_shuff(isnan(X_shuff)) = 0;
    R2_shuff_full = fit_cv_ridge_pooled(X_shuff, y, cv_idx, k_folds, lambdaVal);
    ndr = R2_shuff_full - R2_red;
end

function [pooled_R2_test, y_pred_all] = fit_cv_ridge_pooled(X, y, cv_idx, k_folds, lambdaVal)
    N = length(y); y_pred_all = nan(N, 1); 
    for f = 1:k_folds
        test_idx = (cv_idx == f); train_idx = (cv_idx ~= f & cv_idx > 0); 
        if sum(train_idx) == 0 || sum(test_idx) == 0; continue; end
        X_tr = [ones(sum(train_idx),1), X(train_idx, :)]; X_te = [ones(sum(test_idx),1), X(test_idx, :)];
        penalty = lambdaVal * eye(size(X_tr, 2)); penalty(1, 1) = 0; 
        b = (X_tr' * X_tr + penalty) \ (X_tr' * y(train_idx));
        y_pred_all(test_idx) = X_te * b;
    end
    valid_idx = ~isnan(y_pred_all) & cv_idx > 0;
    if sum(valid_idx) == 0, pooled_R2_test = 0; return; end
    sst = sum((y(valid_idx) - mean(y(valid_idx))).^2); sse = sum((y(valid_idx) - y_pred_all(valid_idx)).^2);
    if sst > 0, pooled_R2_test = 1 - (sse/sst); else, pooled_R2_test = 0; end
end

function [pooled_R2_red, y_pred_all, fold_mse_red, fold_mse_full] = fit_cv_ridge_compare(X_red, X_full, y, cv_idx, k_folds, lambdaVal)
    N = length(y); y_pred_all = nan(N, 1); fold_mse_red = zeros(k_folds, 1); fold_mse_full = zeros(k_folds, 1);
    for f = 1:k_folds
        test_idx = (cv_idx == f); train_idx = (cv_idx ~= f & cv_idx > 0); 
        if sum(train_idx) == 0 || sum(test_idx) == 0; continue; end
        X_tr_red = [ones(sum(train_idx),1), X_red(train_idx, :)]; X_te_red = [ones(sum(test_idx),1), X_red(test_idx, :)];
        pen_red = lambdaVal * eye(size(X_tr_red, 2)); pen_red(1,1) = 0; b_red = (X_tr_red' * X_tr_red + pen_red) \ (X_tr_red' * y(train_idx));
        y_pred_red = X_te_red * b_red; y_pred_all(test_idx) = y_pred_red;
        X_tr_full = [ones(sum(train_idx),1), X_full(train_idx, :)]; X_te_full = [ones(sum(test_idx),1), X_full(test_idx, :)];
        pen_full = lambdaVal * eye(size(X_tr_full, 2)); pen_full(1,1) = 0; b_full = (X_tr_full' * X_tr_full + pen_full) \ (X_tr_full' * y(train_idx));
        fold_mse_red(f) = mean((y(test_idx) - y_pred_red).^2); fold_mse_full(f) = mean((y(test_idx) - X_te_full * b_full).^2);
    end
    valid_idx = ~isnan(y_pred_all) & cv_idx > 0; sst = sum((y(valid_idx) - mean(y(valid_idx))).^2); sse = sum((y(valid_idx) - y_pred_all(valid_idx)).^2);
    if sst > 0, pooled_R2_red = 1 - (sse/sst); else, pooled_R2_red = 0; end
end
