 % 1. Select the session
    [fileName, pathName] = uigetfile('~/Desktop/VIP-GCaMP-sData/*.mat', 'Select a session to inspect');
    if isequal(fileName, 0), return; end
    
    fprintf('Loading %s...\n', fileName);
    S = load(fullfile(pathName, fileName));
    sData = S.sData;
    
    % 2. Extract continuous behavior signals
    time = double(sData.behavior.signals.timeDs(:));
    pos = double(sData.behavior.signals.corridorPositionDs(:));
    
    vel = double(sData.behavior.signals.velocityDs(:));
    % --- PATCH VELOCITY ARTIFACT FIRST --- 
    if length(vel) > 2
        vel(1:2) = vel(3);
    end
    lick = double(sData.behavior.signals.lickIfreqDs(:));
    
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
            pupil = double(pupilDiameterUs(imagingFrameInd));
            pupil = pupil(:);
        else
            pupil = nan(length(time), 1);
        end
    catch ME
        fprintf('Note: Pupil data could not be extracted (%s). Defaulting to NaNs.\n', ME.message);
        pupil = nan(length(time), 1);
    end
    
    % 3. Extract continuous dF/F & Robust Dimension Matching
    dff_matrix = sData.imdata.roiSignals(2).dff;
    
    % Check for accidental transpose
    if size(dff_matrix, 2) ~= length(time) && size(dff_matrix, 1) == length(time)
        dff_matrix = dff_matrix'; 
    end
    
    % Enforce strict dimension matching across all arrays to prevent crashes
    N = min([size(dff_matrix, 2), length(time), length(pos), length(vel), length(lick), length(pupil)]);
    dff_matrix = dff_matrix(:, 1:N);
    time = time(1:N); 
    pos = pos(1:N); 
    vel = vel(1:N); 
    lick = lick(1:N); 
    pupil = pupil(1:N);
    nCells = size(dff_matrix, 1);
    
    % 4. Extract Reward Events
    if isfield(sData, 'daqdata') && isfield(sData.daqdata, 'frameIndex')
        imagingFrameInd = find(diff(sData.daqdata.frameIndex) == 1) + 1;
        N_frames = min(length(imagingFrameInd), length(time));
        imagingFrameInd = imagingFrameInd(1:N_frames);
        
        if isfield(sData.daqdata, 'waterValve')
            valve = double(sData.daqdata.waterValve(imagingFrameInd) > 0.5);
            valve = valve(:); 
        else
            valve = zeros(length(time), 1);
        end
    else
        valve = zeros(length(time), 1);
    end
    
    % 5. Identify Block Transitions
    trialIdx = sData.behavior.signals.trialIndexesDs(1:N);
    
    % Take the most frequent context ID across the entire trial to ignore edge artifacts
    ctx_raw = mode(sData.behavior.trialMatrices.contextMatrix, 2);
    ctx = round(fillmissing(ctx_raw, 'nearest')); 
    
    change_trials = find(diff(ctx) ~= 0) + 1;
    trans_times = [];
    context_ids = [ctx(1)]; 
    
    for i = 1:length(change_trials)
        first_frame = find(trialIdx == change_trials(i), 1, 'first');
        if ~isempty(first_frame)
            t_trans = time(first_frame);
            if t_trans > 10 % Ignore start-up artifacts
                trans_times(end+1) = t_trans;
                context_ids(end+1) = ctx(change_trials(i));
            end
        end
    end
    
    % 6. Setup the Figure
    figure('Name', ['Session Inspector: ', fileName], 'Color', 'w', 'Position', [100, 100, 1000, 1000]);
    
    % --- Panel 1: VIP Cell dF/F Traces (50% OF THE WINDOW) ---
    ax1 = subplot(8, 1, [1 2 3 4]);
    spacing = prctile(dff_matrix(:), 99) * 1.2; 
    if spacing == 0 || isnan(spacing), spacing = 2; end
    
    dff_offset = dff_matrix + (1:nCells)' * spacing;
    plot(time, dff_offset', 'k', 'LineWidth', 0.5);
    
    ylabel('VIP Cells \DeltaF/F');
    title(sprintf('Continuous \\DeltaF/F (%d Cells)', nCells));
    set(gca, 'TickDir', 'out', 'YTick', []); 
    axis tight;
    
    % Add transition lines and clean block labels
    hold on;
    y_max = max(dff_offset(:)) + spacing;
    ylim([0, y_max + spacing*2]); 
    
    block_centers = [];
    if ~isempty(trans_times)
        block_centers(1) = trans_times(1) / 2;
        for i = 1:length(trans_times)-1
            block_centers(i+1) = mean([trans_times(i), trans_times(i+1)]);
        end
        block_centers(end+1) = mean([trans_times(end), max(time)]);
    else
        block_centers(1) = max(time)/2;
    end
    
    for t_idx = 1:length(trans_times)
        xline(trans_times(t_idx), 'k--', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    
    % Label blocks F-N-F or print Context IDs
    for b = 1:length(block_centers)
        if b <= length(context_ids)
            if length(context_ids) == 3 && b == 1, label_txt = 'Familiar';
            elseif length(context_ids) == 3 && b == 2, label_txt = 'Novel Environment';
            elseif length(context_ids) == 3 && b == 3, label_txt = 'Familiar';
            else, label_txt = sprintf('Ctx: %d', context_ids(b));
            end
            
            text(block_centers(b), y_max, label_txt, ...
                'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
        end
    end
    hold off;
    
    % --- Panel 2: VR Position ---
    ax2 = subplot(8, 1, 5);
    plot(time, pos, 'LineWidth', 1.5, 'Color', [0.2 0.6 0.2]);
    ylabel('Position (cm)');
    set(gca, 'TickDir', 'out', 'XTickLabel', []);
    axis tight; ylim([0 260]);
    add_transition_lines(trans_times);
    
    % --- Panel 3: Running Speed ---
    ax3 = subplot(8, 1, 6);
    plot(time, vel, 'LineWidth', 1.5, 'Color', [0.8 0.4 0.2]);
    ylabel('Velocity (cm/s)');
    set(gca, 'TickDir', 'out', 'XTickLabel', []);
    axis tight;
    add_transition_lines(trans_times);
    
    % --- Panel 4: Lick Rate and Rewards ---
    ax4 = subplot(8, 1, 7);
    hold on;
    h_lick = plot(time, lick, 'LineWidth', 1.5, 'Color', [0.2 0.6 0.8], 'DisplayName', 'Lick Rate');
    
    reward_times = time(valve > 0);
    reward_times = reward_times(:); 
    tick_height = max(lick) * 0.15; 
    if isnan(tick_height) || tick_height == 0, tick_height = 1; end
    
    if ~isempty(reward_times)
        reward_onsets = reward_times(diff([0; reward_times]) > 1); 
        for r = 1:length(reward_onsets)
            plot([reward_onsets(r), reward_onsets(r)], [0, tick_height], 'r-', 'LineWidth', 2.5, 'HandleVisibility', 'off');
        end
        h_rew = plot(NaN, NaN, 'r-', 'LineWidth', 2.5, 'DisplayName', 'Reward Delivery');
        legend([h_lick, h_rew], 'Location', 'northeast', 'FontSize', 10);
    else
        legend(h_lick, 'Location', 'northeast', 'FontSize', 10);
    end
    ylabel('Lick Rate');
    set(gca, 'TickDir', 'out', 'XTickLabel', []);
    axis tight;
    add_transition_lines(trans_times);
    hold off;
    
    % --- Panel 5: Pupil Diameter ---
    ax5 = subplot(8, 1, 8);
    hold on;
    if all(isnan(pupil))
        text(mean(time), 0, 'No Pupil Data Found', 'HorizontalAlignment', 'center', 'Color', [0.5 0.5 0.5]);
        set(gca, 'YTick', []);
    else
        plot(time, pupil, 'LineWidth', 1.5, 'Color', [0.4 0.2 0.6]);
    end
    ylabel('Pupil (a.u.)');
    xlabel('Time (seconds)');
    set(gca, 'TickDir', 'out');
    axis tight;
    add_transition_lines(trans_times);
    hold off;
    
    % 7. Link Axes for synchronized zooming
    linkaxes([ax1, ax2, ax3, ax4, ax5], 'x');
    drawnow;


% Helper function (HandleVisibility off prevents messy legends)
function add_transition_lines(trans_times)
    hold on;
    for t_idx = 1:length(trans_times)
        xline(trans_times(t_idx), 'k--', 'LineWidth', 2.5, 'HandleVisibility', 'off');
    end
    hold off;
end