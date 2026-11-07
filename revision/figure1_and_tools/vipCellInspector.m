function vipCellInspector()
    % 1. File Selection
    userHome = getenv('HOME'); 
    folderPath = fullfile(userHome, 'Desktop', 'VIP-GCaMP-sData');
    if ~exist(folderPath, 'dir'), folderPath = fullfile(userHome, 'Desktop'); end
    
    [fileName, pathName] = uigetfile('*.mat', 'Select a VIP session file', folderPath);
    if isequal(fileName, 0), return; end
    
    % Load data
    fprintf('Loading %s...\n', fileName);
    data = load(fullfile(pathName, fileName));
    sData = data.sData;
    
    % Check for VIP neural data
    if isfield(sData.imdata, 'roiSignals') && length(sData.imdata.roiSignals) >= 2
        dff_all = sData.imdata.roiSignals(2).dff;
    else
        error('Could not find VIP neural data in sData.imdata.roiSignals(2).dff');
    end
    nCells = size(dff_all, 1);
    
    % 2. Setup GUI (Slightly taller to fit the dropdown)
    uiFig = uifigure('Name', ['VIP Multi-Feature Inspector: ' fileName], 'Position', [100 100 300 190]);
    
    envDropdown = uidropdown(uiFig, 'Position', [20 130 260 30], ...
        'Items', {'Familiar Context', 'Novel Context'});
        
    spinner = uispinner(uiFig, 'Position', [20 80 260 30], 'Limits', [1 nCells], 'Value', 1);
    
    uibutton(uiFig, 'Text', 'Plot Feature Overlays', 'Position', [20 20 260 40], ...
        'ButtonPushedFcn', @(src, event) plotSpatialFeatures(sData, dff_all, spinner.Value, envDropdown.Value, uiFig));
    
    function plotSpatialFeatures(sData, dff_all, cellID, envChoice, uiFig)
        % 1. Extract Continuous Data & FORCE 1D COLUMN VECTORS
        time = sData.behavior.signals.timeDs(:);
        N = min(size(dff_all, 2), length(time));
        
        dff = zscore(double(dff_all(cellID, 1:N)')); dff = dff(:);
        pos = double(sData.behavior.signals.corridorPositionDs(1:N)); pos = pos(:);
        vel = double(sData.behavior.signals.velocityDs(1:N)); vel = vel(:);
        
        % --- ARTIFACT CATCHER ---
        if length(vel) > 2
            vel(1:2) = vel(3);
        end
        
        acc = gradient(vel); acc = acc(:);
        lick = double(sData.behavior.signals.lickIfreqDs(1:N)); lick = lick(:);
        
        % Robust Valve (Reward) Extraction
        if isfield(sData, 'daqdata') && isfield(sData.daqdata, 'frameIndex')
            img_frames = find(diff(sData.daqdata.frameIndex) == 1) + 1;
            N_frames = min(length(img_frames), N);
            img_frames = img_frames(1:N_frames);
            if isfield(sData.daqdata, 'waterValve')
                valve = double(sData.daqdata.waterValve(img_frames) > 0.5);
                valve = valve(:);
            else
                valve = zeros(N, 1);
            end
        else
            valve = zeros(N, 1);
        end
        if length(valve) < N, valve = [valve; zeros(N - length(valve), 1)]; else, valve = valve(1:N); end
        valve = valve(:);
        
        % --- EXTRACT PUPIL (MATE'S ALIGNMENT METHOD) ---
        try
            if isfield(sData.behavior, 'pupildata') && isfield(sData.behavior.pupildata, 'diameter')
                pupilDiameter = double(sData.behavior.pupildata.diameter);
                pupilDiameter(isnan(pupilDiameter)) = inf; % Preserve NaNs
                
                pupilFrameInd = find(diff(sData.behavior.signals.pupilCamFrameIndexCorr)==1)+1; 
                imagingFrameInd = find(diff(sData.daqdata.frameIndex)==1);
                
                pupilDiameterUs = nan(size(sData.behavior.signals.pupilCamFrameIndexCorr));
                
                % Safely map to the upsampled array
                valid_len = min(length(pupilFrameInd), length(pupilDiameter));
                pupilDiameterUs(pupilFrameInd(1:valid_len)) = pupilDiameter(1:valid_len);
                
                pupilDiameterUs = fillmissing(pupilDiameterUs, 'nearest');
                pupilDiameterUs(isinf(pupilDiameterUs)) = nan;
                
                % Downsample to imaging frames
                pupil_raw = double(pupilDiameterUs(imagingFrameInd));
                
                % Match length N
                if length(pupil_raw) >= N
                    pupil = pupil_raw(1:N);
                else
                    pupil = [pupil_raw(:); nan(N - length(pupil_raw), 1)];
                end
                pupil = pupil(:);
            else
                pupil = nan(N, 1); 
            end
        catch ME
            fprintf('Note: Pupil data could not be extracted (%s). Defaulting to NaNs.\n', ME.message);
            pupil = nan(N, 1);
        end
        
        % 2. Context Masking Logic
        ctx_raw = mode(sData.behavior.trialMatrices.contextMatrix, 2);
        ctx = round(fillmissing(ctx_raw, 'nearest')); ctx = ctx(:);
        trialIdx = sData.behavior.signals.trialIndexesDs(1:N); trialIdx = trialIdx(:);
        
        valid_trials = trialIdx > 0 & trialIdx <= length(ctx);
        continuous_ctx = nan(N, 1);
        continuous_ctx(valid_trials) = ctx(trialIdx(valid_trials));
        continuous_ctx = fillmissing(continuous_ctx, 'nearest');
        
        fam_id = ctx(1);
        unique_ctx = unique(ctx);
        nov_id = unique_ctx(unique_ctx ~= fam_id);
        
        % Apply the selected mask
        if strcmp(envChoice, 'Familiar Context')
            mask = (continuous_ctx == fam_id);
        else
            if isempty(nov_id)
                uialert(uiFig, 'No Novel Context found in this session.', 'Missing Data');
                return;
            end
            mask = (continuous_ctx == nov_id(1));
        end
        mask = mask(:);
        
        % 3. Binning over Position
        numBins = 50;
        binEdges = linspace(0, max(pos), numBins + 1);
        binCenters = linspace(0, max(pos), numBins);
        
        mDFF = nan(numBins,1); sDFF = nan(numBins,1);
        mVel = nan(numBins,1); sVel = nan(numBins,1);
        mAcc = nan(numBins,1); sAcc = nan(numBins,1);
        mLick = nan(numBins,1); sLick = nan(numBins,1);
        mValve = nan(numBins,1);
        mPupil = nan(numBins,1); sPupil = nan(numBins,1);
        
        for b = 1:numBins
            idx = (pos >= binEdges(b) & pos < binEdges(b+1)) & mask;
            n_pts = sum(idx);
            if n_pts > 0
                mDFF(b) = nanmean(dff(idx)); sDFF(b) = nanstd(dff(idx))/sqrt(n_pts);
                mVel(b) = nanmean(vel(idx)); sVel(b) = nanstd(vel(idx))/sqrt(n_pts);
                mAcc(b) = nanmean(acc(idx)); sAcc(b) = nanstd(acc(idx))/sqrt(n_pts);
                mLick(b) = nanmean(lick(idx)); sLick(b) = nanstd(lick(idx))/sqrt(n_pts);
                mValve(b) = nanmean(valve(idx)); 
                mPupil(b) = nanmean(pupil(idx)); sPupil(b) = nanstd(pupil(idx))/sqrt(n_pts);
            end
        end
        
        % 4. Plotting
        f = figure(882); clf(f);
        set(f, 'Color', 'w', 'Name', ['VIP Cell #' num2str(cellID)], 'Position', [200 100 800 850]);
        rewardLines = [100, 220]; 
        
        % Embedded Helper function: Plots shaded errors safely
        function plotShaded(x, y, e, col)
            v = ~isnan(y); 
            if any(v)
                fill([x(v), fliplr(x(v))], [y(v)'+e(v)', fliplr(y(v)'-e(v)')], col, 'EdgeColor', 'none', 'FaceAlpha', 0.3);
            end
        end
        
        % --- Panel 1: DFF ---
        ax1 = subplot(4,1,1); hold on;
        plotShaded(binCenters, mDFF, sDFF, [0.3 0.3 0.3]);
        plot(binCenters, mDFF, 'k', 'LineWidth', 2);
        xline(rewardLines, '--r', 'Alpha', 0.5);
        ylabel('Z-Scored \DeltaF/F'); 
        title(sprintf('VIP Cell #%d (%s)', cellID, envChoice));
        axis tight; set(gca, 'TickDir', 'out', 'XTickLabel', []);
        
        % --- Panel 2: Velocity & Acceleration ---
        ax2 = subplot(4,1,2); hold on;
        yyaxis left 
        plotShaded(binCenters, mVel, sVel, [0.2 0.5 0.9]);
        plot(binCenters, mVel, '-b', 'LineWidth', 2);
        ylabel('Velocity (cm/s)');
        
        yyaxis right 
        plotShaded(binCenters, mAcc, sAcc, [0.9 0.4 0.1]);
        plot(binCenters, mAcc, '-', 'Color', [0.8 0.3 0], 'LineWidth', 1.5);
        ylabel('Acceleration');
        xline(rewardLines, '--k', 'Alpha', 0.3);
        title('Locomotion'); axis tight; set(gca, 'TickDir', 'out', 'XTickLabel', []);
        
        % --- Panel 3: Licking & Valve Opening ---
        ax3 = subplot(4,1,3); hold on;
        yyaxis left 
        plotShaded(binCenters, mLick, sLick, [0.8 0.2 0.8]);
        plot(binCenters, mLick, '-m', 'LineWidth', 2);
        ylabel('Lick Rate (Hz)');
        
        yyaxis right 
        bar(binCenters, mValve, 'FaceColor', [0.2 0.8 0.2], 'EdgeColor', 'none', 'FaceAlpha', 0.2);
        ylabel('Reward Dens.');
        xline(rewardLines, '--k', 'Alpha', 0.3);
        title('Consummatory Behavior'); axis tight; set(gca, 'TickDir', 'out', 'XTickLabel', []);
        
        % --- Panel 4: Pupil Diameter ---
        ax4 = subplot(4,1,4); hold on;
        if all(isnan(mPupil))
            text(mean(binCenters), 0, 'No Pupil Data Found', 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
            set(gca, 'YTick', []);
        else
            plotShaded(binCenters, mPupil, sPupil, [0.2 0.7 0.7]);
            plot(binCenters, mPupil, '-', 'Color', [0 0.6 0.6], 'LineWidth', 2);
            ylabel('Pupil Diameter'); 
        end
        xlabel('Position (cm)');
        xline(rewardLines, '--k', 'Alpha', 0.3);
        title('Arousal (Pupil)'); axis tight; set(gca, 'TickDir', 'out');
        
        linkaxes([ax1, ax2, ax3, ax4], 'x'); 
    end
end