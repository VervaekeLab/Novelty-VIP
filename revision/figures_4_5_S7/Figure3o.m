%% ========================================================================
%% FIGURE 3o: backward shift vs distance to reward
%% ========================================================================
% Left panel : per-group LME spatial trend with binned means +/- SEM.
% Right bars : per-mouse OLS slope and intercept (bars +/- SEM, points = mice),
%              with LME between-group contrasts shown as brackets.
% Stats      : linear mixed-effects model with a random slope per mouse
%              (placeFieldShift ~ rewDistance*GroupCat + (1+rewDistance|mouseID);
%              Satterthwaite df), so contrasts are evaluated across mice.
% ========================================================================
clearvars -except tbl; clc;
rng(42);   % reproducible jitter

alpha_threshold = 0.055;   % significance threshold

%% --- 1. SETUP & DATA CLEANING ---
T = tbl;
groups_display = {'Control', 'ArchT', 'ChrimsonR'};
colors = [0.5 0.5 0.5; 0.0 0.45 0.74; 0.85 0.33 0.1];

T.placeFieldShift = double(T.placeFieldShift);
T.rewDistance = double(T.rewDistance);
raw_g = lower(string(T.group));

T.GroupCat = strings(height(T), 1);
T.GroupCat(contains(raw_g, 'arch')) = "ArchT";
T.GroupCat(contains(raw_g, 'chri')) = "ChrimsonR";
is_ctrl = (contains(raw_g, 'ctrl') | contains(raw_g, 'control')) & (T.GroupCat == "");
T.GroupCat(is_ctrl) = "Control";

T = T(T.GroupCat ~= "", :);
T.GroupCat = categorical(T.GroupCat, groups_display);
T.mouseID = categorical(T.mouseID);

%% --- 2. PER-MOUSE OLS (for the bar plots) ---
mouse_slopes = cell(1,3);
mouse_ints = cell(1,3);
for g = 1:3
    sub_T = T(T.GroupCat == groups_display{g}, :);
    uMice = unique(sub_T.mouseID);
    for m = 1:length(uMice)
        idx = sub_T.mouseID == uMice(m);
        if sum(idx) > 5
            mdl = fitlm(sub_T(idx,:), 'placeFieldShift ~ rewDistance');
            mouse_ints{g}(end+1) = mdl.Coefficients.Estimate(1);
            mouse_slopes{g}(end+1) = mdl.Coefficients.Estimate(2);
        end
    end
end

%% --- 3. LME STATISTICS (random slope per mouse -> across-mice df) ---
lme = fitlme(T, 'placeFieldShift ~ rewDistance * GroupCat + (1 + rewDistance | mouseID)');
FE = lme.Coefficients.Estimate;
names = lme.CoefficientNames;

pS = @(nm) coefTest(lme, double(strcmp(names, nm)), 0, 'DFMethod', 'satterthwaite');
pC = @(nm1, nm2) coefTest(lme, double(strcmp(names, nm1)) - double(strcmp(names, nm2)), 0, 'DFMethod', 'satterthwaite');

idx_int = [find(strcmp(names, '(Intercept)')), find(strcmp(names, 'GroupCat_ArchT')), find(strcmp(names, 'GroupCat_ChrimsonR'))];
idx_slope = [find(strcmp(names, 'rewDistance')), find(strcmp(names, 'rewDistance:GroupCat_ArchT')), find(strcmp(names, 'rewDistance:GroupCat_ChrimsonR'))];

lme_lines_int = [FE(idx_int(1)), FE(idx_int(1))+FE(idx_int(2)), FE(idx_int(1))+FE(idx_int(3))];
lme_lines_slope = [FE(idx_slope(1)), FE(idx_slope(1))+FE(idx_slope(2)), FE(idx_slope(1))+FE(idx_slope(3))];

%% --- 4. FIGURE ---
fig = figure('Units', 'centimeters', 'Position', [2, 2, 12, 9], 'Color', 'w');
tlo = tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');
fs_ax = 9; fs_lab = 10; fs_stars = 14;

get_plot_stars = @(p) repmat('*', 1, (p<=alpha_threshold) + (p<0.01) + (p<0.001));

% --- left panel: LME spatial trends + binned means ---
nexttile(1, [2 2]); hold on;
bin_size = 10;
group_offset = 2.0;
edges = floor(min(T.rewDistance)/10)*10 : bin_size : ceil(max(T.rewDistance)/10)*10;
bin_centers = edges(1:end-1) + bin_size/2;

for g = 1:3
    g_data = T(T.GroupCat == groups_display{g}, :);
    x_range = [min(T.rewDistance) max(T.rewDistance)];
    plot(x_range, lme_lines_int(g) + lme_lines_slope(g) * x_range, 'Color', colors(g,:), 'LineWidth', 2.5);

    m_bin = []; s_bin = []; c_bin = [];
    for b = 1:length(edges)-1
        if b == length(edges)-1
            mask = g_data.rewDistance >= edges(b) & g_data.rewDistance <= edges(b+1);
        else
            mask = g_data.rewDistance >= edges(b) & g_data.rewDistance < edges(b+1);
        end
        if sum(mask) >= 5
            m_bin(end+1) = mean(g_data.placeFieldShift(mask));
            s_bin(end+1) = std(g_data.placeFieldShift(mask)) / sqrt(sum(mask));
            c_bin(end+1) = bin_centers(b);
        end
    end
    current_offset = (g - 2) * group_offset;
    errorbar(c_bin + current_offset, m_bin, s_bin, 'o', ...
        'Color', colors(g,:), 'MarkerFaceColor', colors(g,:), ...
        'LineWidth', 1.5, 'MarkerSize', 5, 'CapSize', 0);
end

yline(0, '--k', 'LineWidth', 1, 'Alpha', 0.5);
xlabel('Distance to Reward (cm)', 'FontSize', fs_lab, 'FontWeight', 'bold');
ylabel('Backward Shift (cm)', 'FontSize', fs_lab, 'FontWeight', 'bold');
title('Spatial Trend', 'FontSize', fs_lab);
set(gca, 'Box', 'off', 'FontSize', fs_ax, 'LineWidth', 1, 'TickDir', 'out');

% --- right bars: per-mouse OLS slope (top) and intercept (bottom), LME contrasts ---
for p = 1:2
    if p == 1
        nexttile(3); hold on;
        vals = mouse_slopes;
        p_arch = pS('rewDistance:GroupCat_ArchT');
        p_chri = pS('rewDistance:GroupCat_ChrimsonR');
        p_arch_chri = pC('rewDistance:GroupCat_ChrimsonR', 'rewDistance:GroupCat_ArchT');
        lbl = 'Slope (cm/cm)'; tit = 'Slope';
    else
        nexttile(6); hold on;
        vals = mouse_ints;
        p_arch = pS('GroupCat_ArchT');
        p_chri = pS('GroupCat_ChrimsonR');
        p_arch_chri = pC('GroupCat_ChrimsonR', 'GroupCat_ArchT');
        lbl = 'Intercept (cm)'; tit = 'Intercept at 0 cm';
    end

    max_data_y = -inf; min_data_y = inf;
    for i = 1:3
        mu = mean(vals{i});
        se = std(vals{i}) / sqrt(length(vals{i}));
        patch([i-0.35 i+0.35 i+0.35 i-0.35], [0 0 mu mu], colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1, 'FaceAlpha', 1);
        plot([i, i], [mu-se, mu+se], '-k', 'LineWidth', 1);
        plot([i-0.1, i+0.1], [mu-se, mu-se], '-k', 'LineWidth', 1);
        plot([i-0.1, i+0.1], [mu+se, mu+se], '-k', 'LineWidth', 1);
        x_jitter = i + (rand(size(vals{i}))-0.5)*0.3;
        scatter(x_jitter, vals{i}, 15, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
        max_data_y = max([max_data_y, mu+se, max(vals{i})]);
        min_data_y = min([min_data_y, mu-se, min(vals{i}), 0]);
    end

    % between-group brackets, stacked
    y_range = max_data_y - min_data_y;
    y_offset = max(y_range * 0.15, 0.02);
    current_level = 1;

    if p_arch <= alpha_threshold
        y_brk = max_data_y + (current_level * 1.5 * y_offset);
        plot([1 1 2 2], [y_brk-y_offset/4, y_brk, y_brk, y_brk-y_offset/4], '-k', 'LineWidth', 1);
        text(1.5, y_brk, get_plot_stars(p_arch), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fs_stars, 'FontWeight', 'bold');
        current_level = current_level + 1;
    end
    if p_arch_chri <= alpha_threshold
        y_brk = max_data_y + (current_level * 1.5 * y_offset);
        plot([2 2 3 3], [y_brk-y_offset/4, y_brk, y_brk, y_brk-y_offset/4], '-k', 'LineWidth', 1);
        text(2.5, y_brk, get_plot_stars(p_arch_chri), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fs_stars, 'FontWeight', 'bold');
        current_level = current_level + 1;
    end
    if p_chri <= alpha_threshold
        y_brk = max_data_y + (current_level * 1.5 * y_offset);
        plot([1 1 3 3], [y_brk-y_offset/4, y_brk, y_brk, y_brk-y_offset/4], '-k', 'LineWidth', 1);
        text(2, y_brk, get_plot_stars(p_chri), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', fs_stars, 'FontWeight', 'bold');
        current_level = current_level + 1;
    end
    if current_level > 1
        current_ylims = ylim;
        ylim([current_ylims(1), max_data_y + (current_level * 1.5 * y_offset) + y_offset]);
    end

    yline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
    xlim([0.3 3.7]);
    ylabel(lbl, 'FontSize', fs_lab, 'FontWeight', 'bold');
    title(tit, 'FontSize', fs_lab);
    if p == 1
        set(gca, 'XColor', 'none');
    else
        xticks(1:3); xticklabels({'Control', 'ArchT', 'ChrimsonR'});
    end
    set(gca, 'FontSize', fs_ax, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');
end
set(findall(fig,'-property','FontName'), 'FontName', 'Arial');

%% --- 5. CONSOLE AUDIT ---
fprintf('\n========================================================\n');
fprintf('  Fig 3o: shift ~ distance * group + (1+distance|mouse)\n');
fprintf('  Satterthwaite df; threshold p <= %.3f\n', alpha_threshold);
fprintf('========================================================\n');
get_console_star = @(p) repmat('*', 1, p <= alpha_threshold);

fprintf('\nSLOPE differences (interaction):\n');
fprintf('  Control vs ArchT      p = %.5f %s\n', pS('rewDistance:GroupCat_ArchT'), get_console_star(pS('rewDistance:GroupCat_ArchT')));
fprintf('  Control vs ChrimsonR  p = %.5f %s\n', pS('rewDistance:GroupCat_ChrimsonR'), get_console_star(pS('rewDistance:GroupCat_ChrimsonR')));
fprintf('  ArchT vs ChrimsonR    p = %.5f %s\n', pC('rewDistance:GroupCat_ChrimsonR','rewDistance:GroupCat_ArchT'), get_console_star(pC('rewDistance:GroupCat_ChrimsonR','rewDistance:GroupCat_ArchT')));

fprintf('\nINTERCEPT differences (group at distance = 0):\n');
fprintf('  Control vs ArchT      p = %.5f %s\n', pS('GroupCat_ArchT'), get_console_star(pS('GroupCat_ArchT')));
fprintf('  Control vs ChrimsonR  p = %.5f %s\n', pS('GroupCat_ChrimsonR'), get_console_star(pS('GroupCat_ChrimsonR')));
fprintf('  ArchT vs ChrimsonR    p = %.5f %s\n', pC('GroupCat_ChrimsonR','GroupCat_ArchT'), get_console_star(pC('GroupCat_ChrimsonR','GroupCat_ArchT')));
fprintf('========================================================\n');
