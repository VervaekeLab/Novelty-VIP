%% ========================================================================
%% FIGURE 3l, m: backward shift of place fields in Novel 1 (BTSP)
%% ========================================================================
% Panel l : mouse-weighted cumulative distribution of backward shift.
% Panel m : suppression of LARGE backward shifts - per-mouse fraction of fields
%           shifting beyond -25 cm, shown for all three groups with brackets.
% Stats   : tail fraction compared across mice on per-mouse values (df ~ mice).
%           Reported: Control vs ArchT, Control vs ChrimsonR, ArchT vs ChrimsonR.
% Note    : ArchT specifically suppresses large backward shifts; the median field
%           is unchanged (see deep-dive). ChrimsonR does not differ significantly.
% ========================================================================
clearvars -except tbl; clc;
rng(0);

%% --- 1. SETUP & DATA CLEANING ---
T = tbl;
groups_display = {'Control', 'ArchT', 'ChrimsonR'};
colors = [0.5 0.5 0.5; 0.0 0.45 0.74; 0.85 0.33 0.1];
bin_size = 2;
thresh_cm = -25;
threshold_bins = thresh_cm / bin_size;

raw_g = lower(string(T.group));
T.GroupCat = strings(height(T), 1);
T.GroupCat(contains(raw_g, 'arch')) = "ArchT";
T.GroupCat(contains(raw_g, 'chri')) = "ChrimsonR";
is_ctrl = (contains(raw_g, 'ctrl') | contains(raw_g, 'control')) & (T.GroupCat == "");
T.GroupCat(is_ctrl) = "Control";
T = T(T.GroupCat ~= "", :);
T.GroupCat = categorical(T.GroupCat, groups_display);
T.mouseID = categorical(T.mouseID);
T.placeFieldShift = double(T.placeFieldShift);
T.IsTail = double(T.placeFieldShift < threshold_bins);

% per-mouse tail fraction (one value per mouse)
mice = unique(T.mouseID);
numMice = length(mice);
mouse_group = repmat(categorical(missing, groups_display), numMice, 1);
mouse_tail_frac = NaN(numMice, 1);
for m = 1:numMice
    idx = T.mouseID == mice(m);
    if any(idx)
        mouse_group(m)     = T.GroupCat(find(idx, 1));
        mouse_tail_frac(m) = mean(T.IsTail(idx)) * 100;
    end
end
valid = ~isnan(mouse_tail_frac);
mouse_group = mouse_group(valid);
mouse_tail_frac = mouse_tail_frac(valid);

%% --- 2. STATISTICS: tail fraction across mice (all three pairwise) ---
Tm = table(mouse_tail_frac, reordercats(removecats(mouse_group),groups_display), ...
           'VariableNames', {'tail','GroupCat'});
lme = fitlme(Tm, 'tail ~ GroupCat');
cn = lme.CoefficientNames; E = fixedEffects(lme); DFM = {'DFMethod','Satterthwaite'};
cr = @(g) double(strcmp(cn, ['GroupCat_' g]));
pairs = {'Control','ArchT'; 'Control','ChrimsonR'; 'ArchT','ChrimsonR'};
pPair = zeros(3,1);
fprintf('<strong>===== Fig 3m: large-shift fraction across mice =====</strong>\n');
for k = 1:3
    h = cr(pairs{k,2}) - cr(pairs{k,1});
    [p,~,~,df2] = coefTest(lme, h, 0, DFM{:});
    pPair(k) = p;
    fprintf('  %-9s vs %-9s : df=%.1f | p=%.4f\n', pairs{k,1}, pairs{k,2}, df2, p);
end

%% --- 3. FIGURE LAYOUT ---
% Left (2 cols): CDF (squeezed).  Right (2 cols): cartoon TOP, bars BELOW.
figure('Units','centimeters','Position',[4 4 15 10],'Color','w','PaperPositionMode','auto');
tlo = tiledlayout(2, 3, 'TileSpacing','compact','Padding','compact');
fs_ax = 9; fs_lab = 10;

%% --- PANEL l: MOUSE-WEIGHTED CUMULATIVE PROBABILITY (both rows, 2 cols) ---
nexttile([2 2]); hold on;
x_cdf = -60:1:40;
for g = 1:3
    g_name = groups_display{g};
    subT = T(T.GroupCat == g_name, :);
    subMice = unique(subT.mouseID);
    mouse_cdfs = zeros(length(subMice), length(x_cdf));
    for m = 1:length(subMice)
        m_data = subT.placeFieldShift(subT.mouseID == subMice(m)) * bin_size;
        mouse_cdfs(m, :) = mean(m_data <= x_cdf, 1);
    end
    mu_cdf = mean(mouse_cdfs, 1);
    sem_cdf = std(mouse_cdfs, 0, 1) / sqrt(length(subMice));
    fill([x_cdf fliplr(x_cdf)], [mu_cdf+sem_cdf fliplr(mu_cdf-sem_cdf)], ...
        colors(g,:), 'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
    plot(x_cdf, mu_cdf, 'Color', colors(g,:), 'LineWidth', 2, 'DisplayName', g_name);
end
xline(0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
xline(thresh_cm, ':k', 'LineWidth', 1.5, 'HandleVisibility', 'off');
xlabel('Backward shift (cm)', 'FontSize', fs_lab, 'FontWeight', 'bold');
ylabel('Cumulative probability', 'FontSize', fs_lab, 'FontWeight', 'bold');
xlim([-60 25]); ylim([0 1]); yticks(0:0.2:1);
legend('Location', 'northwest', 'Box', 'off', 'FontSize', 11);
set(gca, 'FontSize', fs_ax, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');

%% --- TOP RIGHT: TAIL-SUPPRESSION CARTOON (row 1, col 3) ---
nexttile(3); hold on;
x_cart = -5:0.1:5;
y_ctrl = normcdf(x_cart, 0, 1);
y_supp = y_ctrl;
y_supp(x_cart < 0) = normcdf(x_cart(x_cart < 0), 0, 0.5);   % shallower lower tail
plot(x_cart, y_ctrl, 'Color', colors(1,:), 'LineWidth', 2);
plot(x_cart, y_supp, 'Color', colors(2,:), 'LineWidth', 2);
% dashed line marking the tail threshold (a bit further left)
xline(1.0, '--k', 'LineWidth', 1);
% clearer arrow pointing into the lower tail
quiver(-2.0, 0.30, 0.7, -0.18, 0, 'Color', colors(2,:), 'LineWidth', 2.5, ...
       'MaxHeadSize', 2, 'AutoScale','off');
% label: larger, sitting just above the arrow
text(-1.4, 0.62, {'reduction','large backward','shifts'}, 'Color', colors(2,:), ...
     'FontSize', 11, 'FontWeight','bold','FontName','Arial','HorizontalAlignment','center');
xlim([-3 3]); ylim([0 1.05]); axis off;

%% --- BOTTOM RIGHT: THREE-GROUP TAIL-FRACTION BARS + BRACKETS (row 2, col 3) ---
nexttile(6); hold on;
mu = zeros(1,3); se = zeros(1,3); vmax = 0;
for g = 1:3
    v = mouse_tail_frac(mouse_group == groups_display{g});
    mu(g) = mean(v); se(g) = std(v)/sqrt(numel(v));
    vmax = max(vmax, max(v));
    patch([g-0.35 g+0.35 g+0.35 g-0.35], [0 0 mu(g) mu(g)], colors(g,:), ...
          'EdgeColor', 'k', 'LineWidth', 1, 'FaceAlpha', 1);
    plot([g g], [mu(g) mu(g)+se(g)], '-k', 'LineWidth', 1);
    plot([g-0.1 g+0.1], [mu(g)+se(g) mu(g)+se(g)], '-k', 'LineWidth', 1);
    xj = g + (rand(numel(v),1)-0.5)*0.3;
    scatter(xj, v, 22, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);
end

% brackets above the highest data point; Ctrl-ArchT (lower) and Ctrl-Chr (upper)
lvl  = vmax * 1.12;
step = vmax * 0.20;
bk = {[1 2] pPair(1); [1 3] pPair(2)};      % ArchT-Chr removed
for b = 1:2
    pr = bk{b,1}; pv = bk{b,2}; yb = lvl + (b-1)*step;
    plot([pr(1) pr(1) pr(2) pr(2)], [yb-0.02*vmax yb yb yb-0.02*vmax], '-k', 'LineWidth', 1);
    if pv < 0.05
        lab = '*'; if pv<0.01, lab='**'; end; if pv<0.001, lab='***'; end
        % asterisk sits just below the bracket line so it clears the bracket above
        text(mean(pr), yb-0.04*vmax, lab, 'HorizontalAlignment','center','VerticalAlignment','bottom', ...
             'FontSize', 15, 'FontWeight','bold');
    else
        text(mean(pr), yb, 'n.s.', 'HorizontalAlignment','center','VerticalAlignment','bottom', ...
             'FontSize', 9, 'FontName','Arial');
    end
end

xticks(1:3); xticklabels({'Ctrl','ArchT','Chr'});
ylabel('Backward shifts > 25 cm (%)', 'FontSize', fs_lab, 'FontWeight','bold');
ylim([0 lvl + 2.2*step]); xlim([0.4 3.6]);
set(gca, 'FontSize', fs_ax, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off');

set(findall(gcf,'-property','FontName'), 'FontName', 'Arial');
