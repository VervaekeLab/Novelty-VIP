function figure1f()
    % 1. Define folder containing the 10 VIP sessions
    dataDir = '~/Desktop/VIP-GCaMP-sData/';
    files = dir(fullfile(dataDir, '*.mat'));
    
    if isempty(files)
        error('No .mat files found in %s', dataDir);
    end
    
    all_cells_trial_activity = [];
    
    fprintf('Pooling data from %d sessions...\n', length(files));
    
    % 2. Loop through sessions and extract ALL cells
    for f = 1:length(files)
        S = load(fullfile(dataDir, files(f).name));
        sData = S.sData;
        
        dff_matrix = sData.imdata.roiSignals(2).dff;
        trialIdx = sData.behavior.signals.trialIndexesDs;
        
        % Ensure dimensions match [nCells x nFrames]
        if size(dff_matrix, 2) ~= length(trialIdx)
            if size(dff_matrix, 1) == length(trialIdx)
                dff_matrix = dff_matrix'; 
            else
                N = min(size(dff_matrix,2), length(trialIdx));
                dff_matrix = dff_matrix(:, 1:N);
                trialIdx = trialIdx(1:N);
            end
        end
        
        nCells = size(dff_matrix, 1);
        nTrialsInSession = max(trialIdx);
        
        % Find transition to Novel Environment
        ctx = sData.behavior.trialMatrices.contextMatrix(:,1);
        trans_idx = find(diff(ctx) ~= 0);
        
        if isempty(trans_idx); continue; end
        trans1 = trans_idx(1) + 1; 
        
        % Calculate mean dF/F per trial
        mean_trial_dff = nan(nCells, nTrialsInSession);
        for t = 1:nTrialsInSession
            frames = (trialIdx == t);
            if any(frames)
                mean_trial_dff(:, t) = mean(dff_matrix(:, frames), 2, 'omitnan');
            end
        end
        
        % Align to 110 trial window (30 Fam, 60 Nov, 20 Fam)
        start_trial = trans1 - 30;
        aligned_dff = nan(nCells, 110);
        
        for t_col = 1:110
            actual_t = start_trial + t_col - 1;
            if actual_t >= 1 && actual_t <= nTrialsInSession
                aligned_dff(:, t_col) = mean_trial_dff(:, actual_t);
            end
        end
        
        % Visually patch NaNs at the edges of short sessions
        aligned_dff = fillmissing(aligned_dff, 'nearest', 2);
        all_cells_trial_activity = [all_cells_trial_activity; aligned_dff];
    end
    
    nTotalCells = size(all_cells_trial_activity, 1);
    fprintf('\nSuccessfully pooled %d VIP cells total.\n', nTotalCells);
    if nTotalCells == 0, error('No cells were extracted.'); end
    
    % 3. Statistics and Sorting (Novelty Response)
    p_values = ones(nTotalCells, 1);
    mean_diff = zeros(nTotalCells, 1);
    
    for c = 1:nTotalCells
        fam_baseline = all_cells_trial_activity(c, 21:30);
        nov_response = all_cells_trial_activity(c, 31:40);
        
        valid = ~isnan(fam_baseline) & ~isnan(nov_response);
        if sum(valid) >= 3
            try
                % Right-tailed t-test for significant upregulation
                [~, p_values(c)] = ttest(nov_response(valid), fam_baseline(valid), 'Tail', 'right');
            catch
                p_values(c) = 1;
            end
        end
        mean_diff(c) = mean(nov_response(valid)) - mean(fam_baseline(valid));
    end
    
    % Sort for Heatmap
    [~, sort_idx] = sort(mean_diff, 'descend');
    sorted_population_dff = all_cells_trial_activity(sort_idx, :);
    
    % Calculate proportions for pie chart
    num_pos = sum(p_values < 0.05);
    pct_pos = round((num_pos / nTotalCells) * 100);
    pct_ns = 100 - pct_pos;
    
    % 4. Calculate Z-scored Population Mean & SEM (For Figure 1g)
    % Z-score each cell's entire 110-trial trace
    cell_means = mean(all_cells_trial_activity, 2, 'omitnan');
    cell_stds = std(all_cells_trial_activity, 0, 2, 'omitnan');
    z_scored_traces = (all_cells_trial_activity - cell_means) ./ cell_stds;
    
    pop_mean = mean(z_scored_traces, 1, 'omitnan');
    pop_sem = std(z_scored_traces, 0, 1, 'omitnan') ./ sqrt(sum(~isnan(z_scored_traces), 1));
    
    % 5. Fit Exponential Decay to the first 45 novel trials (31 to 75)
    x_fit = (1:45)';
    y_data = pop_mean(31:75)';
    % Objective: y = A * exp(-x/tau) + C
    obj_fun = @(p) sum((p(1) * exp(-x_fit / p(2)) + p(3) - y_data).^2);
    % Initial guesses: [Amplitude, tau, Baseline]
    p0 = [max(y_data)-min(y_data), 10, min(y_data)];
    p_opt = fminsearch(obj_fun, p0, optimset('Display','off'));
    
    fit_curve = p_opt(1) * exp(-x_fit / p_opt(2)) + p_opt(3);
    tau_val = p_opt(2);
    
    % ---------------------------------------------------------
    % 6. Plotting
    % ---------------------------------------------------------
    figure('Name', 'Figure 1f and 1g', 'Color', 'w', 'Position', [100, 100, 600, 700]);
    
    % --- TOP PANEL: Heatmap (Figure 1f) ---
    subplot(3, 1, [1 2]); 
    ax1 = gca;
    set(ax1, 'Color', [0.9 0.9 0.9]); % Grey background for NaNs
    
    im = imagesc(sorted_population_dff);
    im.AlphaData = ~isnan(sorted_population_dff);
    colormap(ax1, 'parula'); 
    
    c_bar = colorbar;
    c_bar.Label.String = '\DeltaF/F';
    c_bar.Label.FontSize = 12;
    clim([0, 2]); 
    
    set(gca, 'TickDir', 'out', 'FontSize', 12);
    ylabel('Cells');
    xticks(0:10:110); xticklabels([]); % Hide X labels, share with bottom
    yticks([1, 20:20:nTotalCells]);
    
    hold on;
    xline(30.5, 'w--', 'LineWidth', 2);
    xline(90.5, 'w--', 'LineWidth', 2);
    
    % Top Labels
    text(15, -4, 'Familiar', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Clipping', 'off');
    text(60, -4, 'Novel environment', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Clipping', 'off');
    text(100, -4, 'Familiar', 'HorizontalAlignment', 'center', 'FontSize', 14, 'Clipping', 'off');
    hold off;
    
    % --- BOTTOM PANEL: Mean Trace (Figure 1g) ---
    subplot(3, 1, 3);
    hold on;
    
    % Shaded Background for Novel Environment
    fill([30.5 90.5 90.5 30.5], [-2 -2 3 3], [0.92 0.92 0.92], 'EdgeColor', 'none');
    text(33, 1.3, 'Novel env.', 'FontSize', 14);
    
    % Plot SEM Shading
    trials_x = 1:110;
    fill([trials_x, fliplr(trials_x)], [pop_mean + pop_sem, fliplr(pop_mean - pop_sem)], ...
         [0.7 0.7 0.7], 'EdgeColor', 'none', 'FaceAlpha', 0.6);
     
    % Plot Population Mean
    plot(trials_x, pop_mean, 'k-', 'LineWidth', 2.5);
    
    % Plot Exponential Fit
    plot(31:75, fit_curve, 'r--', 'LineWidth', 2.5);
    text(35, fit_curve(end) - 0.2, sprintf('\\tau = %.1f trials', tau_val), 'Color', [0.8 0.1 0.1], 'FontSize', 14, 'FontWeight', 'bold');
    
    % Formatting
    set(gca, 'TickDir', 'out', 'FontSize', 12, 'Layer', 'top');
    xlim([0 110]);
    ylim([-0.6 1.6]);
    xticks(0:20:100);
    xlabel('Trial number');
    ylabel('Z-scored \DeltaF/F');
    hold off;
    
    % --- INSET: Pie Chart ---
    % Creating a floating axes for the pie chart in the top right of the bottom panel
    ax_pie = axes('Position', [0.73 0.22 0.15 0.15]); 
    p = pie(ax_pie, [pct_ns, pct_pos]);
    
    % Style pie chart
    p(1).FaceColor = [0.3 0.3 0.3]; % Dark grey
    p(3).FaceColor = [0.95 0.7 0.4]; % Light orange
    p(1).EdgeColor = 'k'; p(3).EdgeColor = 'k';
    
    % Custom text labels inside the pie
    p(2).String = sprintf('%d%%\nn.s.', pct_ns);
    p(2).Color = 'w'; p(2).FontSize = 10; p(2).FontWeight = 'bold';
    
    p(4).String = sprintf('%d%%\npos.', pct_pos);
    p(4).Color = 'w'; p(4).FontSize = 10; p(4).FontWeight = 'bold';
    
    title(ax_pie, 'Modulation', 'FontSize', 12);
end