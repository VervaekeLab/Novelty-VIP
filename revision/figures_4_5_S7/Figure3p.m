%% ========================================================================
%% FIGURE 3p: place-field width vs running speed
%% ========================================================================
% Main panel : per-group width-vs-speed trend (per-mouse-binned means +/- SEM)
%              with OLS trend lines.
% Inset      : per-mouse width-vs-speed slope (bars +/- SEM, points = mice).
% Stats      : linear mixed-effects model with a random slope per mouse
%              (widthData ~ speedC*GroupCat + (1+speedC|mouseID); Satterthwaite df),
%              so contrasts are evaluated across mice. Speed is mean-centred so the
%              intercept contrasts compare width at mean speed (not extrapolated).
%              Inset within-group asterisks mark slope ~= 0; between-group slope
%              differences are shown as brackets.
% ========================================================================
clearvars -except tbl; clc;
rng(42);   % reproducible jitter

%% --- 1. SETUP & DATA PREP ---
bin_size = 10;
groups_display = {'Control', 'ArchT', 'ChrimsonR'};
colors = [0.1 0.1 0.1; 0.0 0.45 0.74; 0.85 0.33 0.1];
group_offset = 1.0;

T = tbl;
T.widthData = double(T.widthData);
T.speedData = double(T.speedData);

raw_g = lower(string(T.group));
T.GroupCat = strings(height(T), 1);
T.GroupCat(contains(raw_g, 'arch')) = "ArchT";
T.GroupCat(contains(raw_g, 'chri')) = "ChrimsonR";
is_ctrl = (contains(raw_g, 'ctrl') | contains(raw_g, 'control')) & (T.GroupCat == "");
T.GroupCat(is_ctrl) = "Control";
T = T(T.GroupCat ~= "", :);
T.GroupCat = categorical(T.GroupCat, groups_display);
T.GroupCat = reordercats(T.GroupCat, groups_display);
T.mouseID = categorical(T.mouseID);

%% --- 2. LME STATISTICS (random slope per mouse -> across-mice df) ---
T.speedC = T.speedData - mean(T.speedData, 'omitnan');   % centre so intercept = width at mean speed
lme = fitlme(T, 'widthData ~ speedC * GroupCat + (1 + speedC | mouseID)');
names = lme.CoefficientNames;
cr = @(tok) coefRow(names, tok);   % order-independent coefficient selector (local fn)

% Random-slope SD (if ~0, slope p-values are unreliable)
psi = covarianceParameters(lme); sdSlope = NaN;
for ii = 1:numel(psi)
    M = psi{ii};
    if isnumeric(M) && isequal(size(M),[2 2]), sdSlope = sqrt(M(2,2)); end
end

% Within-group slope vs 0
H_ctrl = cr({'speedC'});
H_arch = H_ctrl + cr({'speedC','GroupCat_ArchT'});
H_chri = H_ctrl + cr({'speedC','GroupCat_ChrimsonR'});
p_abs = [coefTest(lme, H_ctrl, 0, 'DFMethod','Satterthwaite'), ...
         coefTest(lme, H_arch, 0, 'DFMethod','Satterthwaite'), ...
         coefTest(lme, H_chri, 0, 'DFMethod','Satterthwaite')];

pS = @(H) coefTest(lme, H, 0, 'DFMethod','Satterthwaite');
% Between-group slope differences (interaction terms)
p_slope = struct();
p_slope.CvA = pS(cr({'speedC','GroupCat_ArchT'}));
p_slope.CvC = pS(cr({'speedC','GroupCat_ChrimsonR'}));
p_slope.AvC = pS(cr({'speedC','GroupCat_ChrimsonR'}) - cr({'speedC','GroupCat_ArchT'}));
% Between-group intercept differences (width at mean speed)
p_int = struct();
p_int.CvA = pS(cr({'GroupCat_ArchT'}));
p_int.CvC = pS(cr({'GroupCat_ChrimsonR'}));
p_int.AvC = pS(cr({'GroupCat_ChrimsonR'}) - cr({'GroupCat_ArchT'}));

%% --- 3. PER-MOUSE OLS SLOPES (for the inset) ---
mouse_slopes = cell(1,3);
for g = 1:3
    sub_T = T(T.GroupCat == groups_display{g}, :);
    uMice = unique(sub_T.mouseID);
    for m = 1:length(uMice)
        idx = sub_T.mouseID == uMice(m);
        if sum(idx) > 5
            mdl = fitlm(sub_T(idx,:), 'widthData ~ speedData');
            mouse_slopes{g}(end+1) = mdl.Coefficients.Estimate(2);
        end
    end
end

%% --- 4. FIGURE: MAIN PANEL + INSET ---
fig_width = 12; fig_height = 9;
figure('Units', 'centimeters', 'Position', [10, 10, fig_width, fig_height], 'Color', 'w');
ax_main = axes('Position', [0.15 0.15 0.8 0.8]); hold(ax_main, 'on');

edges = 5 : bin_size : 75;
bin_centers = edges(1:end-1) + bin_size/2;

for g = 1:3
    grp_data = T(T.GroupCat == groups_display{g}, :);
    if isempty(grp_data), continue; end
    uMice = unique(grp_data.mouseID);

    plot_means = []; plot_sems = []; plot_x = [];
    for b = 1:length(edges)-1
        m_means = [];
        for m = 1:length(uMice)
            mask = (grp_data.mouseID == uMice(m)) & (grp_data.speedData >= edges(b)) & (grp_data.speedData < edges(b+1));
            if sum(mask) >= 1, m_means(end+1) = mean(grp_data.widthData(mask)); end
        end
        n_mice = length(m_means);
        if n_mice >= 2
            plot_means(end+1) = mean(m_means);
            plot_sems(end+1)  = std(m_means) / sqrt(n_mice);
            plot_x(end+1)     = bin_centers(b);
        end
    end

    errorbar(ax_main, plot_x + (g-2)*group_offset, plot_means, plot_sems, 'o-', ...
        'Color', colors(g,:), 'MarkerFaceColor', 'w', 'MarkerEdgeColor', colors(g,:), ...
        'LineWidth', 1.5, 'MarkerSize', 5, 'CapSize', 0, 'DisplayName', groups_display{g});

    mdl = fitlm(grp_data, 'widthData ~ speedData');
    beta = mdl.Coefficients.Estimate;
    plot(ax_main, [min(edges), max(edges)], beta(1) + beta(2) * [min(edges), max(edges)], '--', 'Color', colors(g,:), 'LineWidth', 1.5, 'HandleVisibility', 'off');
end

legend(ax_main, 'Location', 'northwest', 'Box', 'off', 'FontSize', 9);
xlabel(ax_main, 'Speed during induction (cm/s)', 'FontSize', 11, 'FontWeight', 'bold');
ylabel(ax_main, 'Place field width (cm)', 'FontSize', 11, 'FontWeight', 'bold');
xlim(ax_main, [5 75]); xticks(ax_main, 10:10:70); ylim(ax_main, [0 65]);
set(ax_main, 'FontSize', 10, 'LineWidth', 1.5, 'TickDir', 'out', 'Box', 'off');

% --- inset: per-mouse slopes ---
ax_inset = axes('Position', [0.60 0.17 0.35 0.42]); hold(ax_inset, 'on');

data_top = -inf;
for i = 1:3
    mu = mean(mouse_slopes{i}); se = std(mouse_slopes{i}) / sqrt(length(mouse_slopes{i}));
    data_top = max([data_top, mu+se, max(mouse_slopes{i})]);
end
star_y_level = data_top + 0.15;   % uniform within-group asterisk height

for i = 1:3
    mu = mean(mouse_slopes{i}); se = std(mouse_slopes{i}) / sqrt(length(mouse_slopes{i}));
    patch(ax_inset, [i-0.35 i+0.35 i+0.35 i-0.35], [0 0 mu mu], colors(i,:), 'EdgeColor', 'k', 'LineWidth', 1, 'FaceAlpha', 1);
    plot(ax_inset, [i, i], [mu-se, mu+se], '-k', 'LineWidth', 1);
    plot(ax_inset, [i-0.1, i+0.1], [mu-se, mu-se], '-k', 'LineWidth', 1);
    plot(ax_inset, [i-0.1, i+0.1], [mu+se, mu+se], '-k', 'LineWidth', 1);
    scatter(ax_inset, i + (rand(size(mouse_slopes{i}))-0.5)*0.3, mouse_slopes{i}, 20, 'w', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

    if p_abs(i) < 0.05   % within-group slope ~= 0
        text(ax_inset, i, star_y_level, '*', 'FontSize', 15, 'FontWeight', 'bold', 'HorizontalAlignment', 'center', 'Color', colors(i,:));
    end
end

% between-group slope brackets
bracket_y_level = star_y_level + 0.25;
bpairs = {[1 2], p_slope.CvA; [2 3], p_slope.AvC; [1 3], p_slope.CvC};
for k = 1:3
    pr = bpairs{k,1}; pv = bpairs{k,2};
    if pv < 0.05
        plot(ax_inset, [pr(1) pr(1) pr(2) pr(2)], [bracket_y_level-0.05, bracket_y_level, bracket_y_level, bracket_y_level-0.05], '-k', 'LineWidth', 1);
        text(ax_inset, mean(pr), bracket_y_level, '*', 'FontSize', 14, 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
    end
end

yline(ax_inset, 0, '--k', 'LineWidth', 1, 'HandleVisibility', 'off');
xlim(ax_inset, [0.3 3.7]); set(ax_inset, 'XColor', 'none');
ylabel(ax_inset, 'Slope (s)', 'FontSize', 11, 'FontWeight', 'bold');
ylim(ax_inset, [-0.4, bracket_y_level + 0.4]);
set(ax_inset, 'FontSize', 10, 'LineWidth', 1, 'TickDir', 'out', 'Box', 'off', 'Color', 'w');

set(findall(gcf,'-property','FontName'), 'FontName', 'Arial');

%% --- 5. CONSOLE AUDIT ---
sig = @(p) ternary(p<0.05,'  <-- significant','');
fprintf('\n=================================================================\n');
fprintf('  Fig 3p: width ~ speed * group + (1+speed|mouse), Satterthwaite\n');
fprintf('  random-slope SD = %.4g%s\n', sdSlope, ternary(sdSlope<1e-3,'  (NEAR-SINGULAR: slope p-values unreliable)',''));
fprintf('=================================================================\n');
fprintf('\nWithin-group slope ~= 0:\n');
for i = 1:3, fprintf('  %-10s p = %.5f%s\n', groups_display{i}, p_abs(i), sig(p_abs(i))); end
fprintf('\nBetween-group SLOPE differences:\n');
fprintf('  Control vs ArchT      p = %.5f%s\n', p_slope.CvA, sig(p_slope.CvA));
fprintf('  Control vs ChrimsonR  p = %.5f%s\n', p_slope.CvC, sig(p_slope.CvC));
fprintf('  ArchT   vs ChrimsonR  p = %.5f%s\n', p_slope.AvC, sig(p_slope.AvC));
fprintf('\nBetween-group INTERCEPT differences (width at mean speed):\n');
fprintf('  Control vs ArchT      p = %.5f%s\n', p_int.CvA, sig(p_int.CvA));
fprintf('  Control vs ChrimsonR  p = %.5f%s\n', p_int.CvC, sig(p_int.CvC));
fprintf('  ArchT   vs ChrimsonR  p = %.5f%s\n', p_int.AvC, sig(p_int.AvC));
fprintf('=================================================================\n');

%% ===================== LOCAL FUNCTIONS =====================
function o = ternary(c,a,b), if c, o=a; else, o=b; end, end

function v = coefRow(names, tokens)
% Order-independent selector: 1 where a coefficient's colon-split parts equal
% the requested token set (handles 'a:b' vs 'b:a' naming), else 0.
    tokens = string(tokens); v = zeros(1, numel(names));
    for k = 1:numel(names)
        parts = string(strsplit(names{k}, ':'));
        if numel(parts)==numel(tokens) && isempty(setdiff(parts,tokens)) ...
                                       && isempty(setdiff(tokens,parts))
            v(k) = 1; return;
        end
    end
    warning('coefRow:noMatch', 'No coefficient matched: %s', strjoin(cellstr(tokens),' & '));
end
