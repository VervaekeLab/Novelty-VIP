%% ========================================================================
%% FIGURE 3n: place-field density vs distance to reward
%% ========================================================================
% Main panel : per-group density-vs-distance trend with binned means +/- SEM.
% Inset      : per-mouse density-vs-distance slope (bars +/- SEM, points = mice).
% Stats      : linear mixed-effects model with a random slope per mouse
%              (Density ~ Distance*Group + (1+Distance|MouseID); Satterthwaite df),
%              so contrasts are evaluated across mice. Inset asterisks mark a
%              within-group gradient (slope ~= 0); between-group differences are
%           d   printed in the console audit.
% ========================================================================
clearvars -except tbl; clc;
rng(0);   % reproducible jitter

%% --- 1. SETUP & DATA CLEANING ---
T = tbl;
binSize = 10;
groups_display = {'Control', 'ArchT', 'ChrimsonR'};
colors = [0.5 0.5 0.5; 0.0 0.45 0.74; 0.85 0.33 0.1];
group_offset = 3.0;   % x-offset for plotting overlapping error bars

raw_g = lower(string(T.group));
T.GroupCat = strings(height(T), 1);
T.GroupCat(contains(raw_g, 'arch')) = "ArchT";
T.GroupCat(contains(raw_g, 'chri')) = "ChrimsonR";
is_ctrl = (contains(raw_g, 'ctrl') | contains(raw_g, 'control')) & (T.GroupCat == "");
T.GroupCat(is_ctrl) = "Control";

T = T(T.GroupCat ~= "", :);
T.GroupCat = categorical(T.GroupCat, groups_display);
T.GroupCat = reordercats(T.GroupCat, groups_display);

edges = floor(min(T.rewDistance)/binSize)*binSize : binSize : ceil(max(T.rewDistance)/binSize)*binSize;
centers = (edges(1:end-1) + binSize/2)';

%% --- 2. PER-MOUSE DENSITIES & SLOPES ---
LME_Table = table();
mouse_slopes = cell(1,3);
all_densities_by_group = cell(1,3);

for g = 1:3
    g_name = groups_display{g};
    sub_T = T(T.GroupCat == g_name, :);
    uMice = unique(sub_T.mouseID);

    g_slopes = [];
    g_densities = [];

    for m = 1:length(uMice)
        dists = sub_T.rewDistance(sub_T.mouseID == uMice(m));
        if length(dists) > 5
            counts = histcounts(dists, edges);
            density = (counts / sum(counts)) * 100;

            tmp = table(repmat(uMice(m), length(centers), 1), ...
                        repmat(string(g_name), length(centers), 1), ...
                        centers, density', ...
                        'VariableNames', {'MouseID', 'Group', 'Distance', 'Density'});
            LME_Table = [LME_Table; tmp]; %#ok<AGROW>

            g_densities(end+1, :) = density; %#ok<AGROW>
            p = polyfit(centers, density', 1);
            g_slopes = [g_slopes; p(1)]; %#ok<AGROW>
        end
    end
    all_densities_by_group{g} = g_densities;
    mouse_slopes{g} = g_slopes;
end

%% --- 3. LME STATISTICS (random slope per mouse -> across-mice df) ---
LME_Table.Group = categorical(LME_Table.Group, groups_display);
LME_Table.Group = reordercats(LME_Table.Group, groups_display);
LME_Table.MouseID = categorical(LME_Table.MouseID);

lme = fitlme(LME_Table, 'Density ~ Distance * Group + (1 + Distance | MouseID)');
names = lme.CoefficientNames;

idx_dist     = strcmp(names, 'Distance');
idx_arch_int = contains(names, 'Distance') & contains(names, 'ArchT');
idx_chri_int = contains(names, 'Distance') & contains(names, 'ChrimsonR');

% Contrast vectors for each group's absolute slope
H_ctrl = double(idx_dist);
H_arch = H_ctrl + double(idx_arch_int);
H_chri = H_ctrl + double(idx_chri_int);

pS_zero = @(H) coefTest(lme, H, 0, 'DFMethod', 'satterthwaite');

% Within-group gradient (slope ~= 0)
p_vs_zero = [pS_zero(H_ctrl), pS_zero(H_arch), pS_zero(H_chri)];

%% --- 4. FIGURE: MAIN PANEL + INSET ---
fig_width = 12; fig_height = 9;
figure('Units', 'centimeters', 'Position', [10, 10, fig_width, fig_height], 'Color', 'w');
ax_main = axes('Position', [0.15 0.15 0.8 0.8]); hold(ax_main, 'on');

% --- main panel: group trend lines + binned means ---
group_maxes = zeros(1,3);
for g = 1:3
    g_densities = all_densities_by_group{g};
    n_mice_g = size(g_densities, 1);

    plot(ax_main, centers, polyval(polyfit(centers, mean(g_densities,1)', 1), centers), ...
        '-', 'Color', colors(g,:), 'LineWidth', 2.5, 'DisplayName', groups_display{g});

    current_offset = (g - 2) * group_offset;
    raw_mean = mean(g_densities, 1);
    raw_sem  = std(g_densities, 0, 1) / sqrt(n_mice_g);
    group_maxes(g) = max(raw_mean + raw_sem);

    errorbar(ax_main, centers + current_offset, raw_mean, raw_sem, 'o', ...
        'Color', colors(g,:), 'MarkerFaceColor', colors(g,:), ...
        'MarkerSize', 5, 'CapSize', 0, 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

xline(ax_main, 0, '--k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel(ax_main, 'Distance to Reward (cm)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ax_main, 'Place Field Density (%)', 'FontSize', 11, 'FontWeight', 'bold');
ylim(ax_main, [0 max(group_maxes)*1.4]);
xlim(ax_main, [min(centers)-12 12]);
set(ax_main, 'FontSize', 10, 'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off');

% --- inset: per-mouse slopes ---
ax_inset = axes('Position', [0.28 0.68 0.24 0.25]); hold(ax_inset, 'on');

global_max_y = -inf;
for i = 1:3
    mu = mean(mouse_slopes{i});
    se = std(mouse_slopes{i}) / sqrt(length(mouse_slopes{i}));
    global_max_y = max([global_max_y, mu+se, max(mouse_slopes{i})]);
end
star_level = global_max_y + 0.02;   % uniform asterisk height

for i = 1:3
    x_pos = i;
    y_vals = mouse_slopes{i};
    mu = mean(y_vals);
    se = std(y_vals) / sqrt(length(y_vals));

    patch(ax_inset, [x_pos-0.35 x_pos+0.35 x_pos+0.35 x_pos-0.35], [0 0 mu mu], colors(i,:), ...
          'EdgeColor', 'k', 'LineWidth', 1, 'FaceAlpha', 1);
    plot(ax_inset, [x_pos, x_pos], [mu-se, mu+se], '-k', 'LineWidth', 1);
    plot(ax_inset, [x_pos-0.1, x_pos+0.1], [mu-se, mu-se], '-k', 'LineWidth', 1);
    plot(ax_inset, [x_pos-0.1, x_pos+0.1], [mu+se, mu+se], '-k', 'LineWidth', 1);
    x_jitter = x_pos + (rand(size(y_vals))-0.5)*0.3;
    scatter(ax_inset, x_jitter, y_vals, 20, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

    if p_vs_zero(i) <= 0.05   % within-group gradient
        text(ax_inset, x_pos, star_level, '*', ...
            'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', colors(i,:));
    end
end

yline(ax_inset, 0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
xticks(ax_inset, 1:3);
xticklabels(ax_inset, {'Ctrl', 'ArchT', 'Chrim'});
ylabel(ax_inset, 'Slope (%/cm)', 'FontSize', 11, 'FontWeight', 'bold');
current_ylims = ylim(ax_inset);
ylim(ax_inset, [current_ylims(1), star_level + 0.02]);
set(ax_inset, 'FontSize', 10, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off', 'Color', 'w');

set(findall(gcf,'-property','FontName'), 'FontName', 'Arial');

%% --- 5. CONSOLE AUDIT ---
fprintf('\n================================================================\n');
fprintf('  Fig 3n: recruitment slopes (LME, Satterthwaite)\n');
fprintf('================================================================\n');
fprintf('\nWithin-group gradient (slope ~= 0):\n');
for i = 1:3
    fprintf('  %-10s : p = %.5f\n', groups_display{i}, p_vs_zero(i));
end
fprintf('\nBetween-group differences (interaction):\n');
pS_diff = @(H_vec) coefTest(lme, H_vec, 0, 'DFMethod', 'satterthwaite');
fprintf('  ArchT     vs Control : p = %.5f\n', pS_diff(double(idx_arch_int)));
fprintf('  ChrimsonR vs Control : p = %.5f\n', pS_diff(double(idx_chri_int)));
fprintf('================================================================\n');
