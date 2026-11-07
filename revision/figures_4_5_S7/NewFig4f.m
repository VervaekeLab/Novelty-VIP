%% FIGURE 4e/f - PLACE-FIELD DRIFT: GROUP COMPARISON ACROSS MICE
% Data: tbl in workspace with shiftData (drift, cm), group, mouseID.
%
% Design: group is a between-mouse factor, so drift is averaged within each
% mouse and compared across mice with a linear model on the per-mouse values
% (fitlm, meanDrift ~ group, pairwise contrasts). Residual df = n_mice - 3.
% All tests two-sided.
%
% Why not fitlme: with no random effects, fitlme fits by maximum likelihood
% and estimates the residual variance as SSE/n instead of SSE/(n-p). With
% Satterthwaite this returns df = n (15) and p-values that are too small.
% fitlm gives the exact test. The old fit is printed below for comparison.

clearvars -except tbl; clc;

%% 1. Setup
T = tbl;
gs = lower(string(T.group)); G = strings(height(T),1);
G(contains(gs,'ctrl')) = "Control";
G(contains(gs,'arch')) = "ArchT";
G(contains(gs,'chri')) = "ChrimsonR";
assert(all(G ~= ""), 'Unrecognised group label(s) in tbl.group.');
grpOrder   = {'Control','ArchT','ChrimsonR'};
T.GroupCat = categorical(G, grpOrder);
T.mouseID  = removecats(categorical(T.mouseID));

%% 2. Per-mouse summaries (compare these with the dots in Fig. 4f)
mice = categories(T.mouseID); nM = numel(mice);
pmMean = nan(nM,1); pmMed = nan(nM,1); pmN = zeros(nM,1); pmGrp = strings(nM,1);
for i = 1:nM
    s = T(T.mouseID == mice{i}, :);
    g = unique(s.GroupCat);
    assert(numel(g) == 1, 'Mouse %s appears in more than one group.', mice{i});
    pmGrp(i)  = string(g);
    pmMean(i) = mean(s.shiftData, 'omitnan');
    pmMed(i)  = median(s.shiftData, 'omitnan');
    pmN(i)    = sum(~isnan(s.shiftData));
end
Tm = table(string(mice), categorical(pmGrp, grpOrder), pmN, pmMean, pmMed, ...
    'VariableNames', {'mouse','GroupCat','nFields','meanDrift','medianDrift'});
Tm = sortrows(Tm, {'GroupCat','meanDrift'});

fprintf('Per-mouse values (compare with the dots in Fig. 4f):\n');
disp(Tm);
fprintf('Group summary across mice:\n');
for g = 1:3
    x = Tm.meanDrift(Tm.GroupCat == grpOrder{g});
    fprintf('  %-9s n = %d | mean = %5.2f | SD = %5.2f | range %5.2f to %5.2f cm\n', ...
        grpOrder{g}, numel(x), mean(x), std(x), min(x), max(x));
end

%% 3. Primary test and robustness battery
pairs = {'Control','ArchT'; 'Control','ChrimsonR'; 'ArchT','ChrimsonR'};
nP = size(pairs,1);

mdl = fitlm(Tm, 'meanDrift ~ GroupCat');            % exact, df = n_mice - 3
[pPri, dPri] = pairContrasts(mdl, pairs);
pHolm = holm(pPri);

mdlMed = fitlm(Tm, 'medianDrift ~ GroupCat');       % per-mouse medians
pMed = pairContrasts(mdlMed, pairs);

pT = nan(nP,1); pW = pT; pRS = pT; pLOO = pT; looMouse = strings(nP,1);
for k = 1:nP
    x = Tm.meanDrift(Tm.GroupCat == pairs{k,1});
    y = Tm.meanDrift(Tm.GroupCat == pairs{k,2});
    [~, pT(k)] = ttest2(x, y);                        % two groups only, pooled
    [~, pW(k)] = ttest2(x, y, 'Vartype', 'unequal');  % Welch
    pRS(k)     = ranksum(x, y);                       % exact at these n
    % leave-one-mouse-out: drop each mouse of the pair, refit primary model
    idx = find(Tm.GroupCat == pairs{k,1} | Tm.GroupCat == pairs{k,2});
    for j = idx'
        m  = fitlm(Tm(setdiff(1:height(Tm), j), :), 'meanDrift ~ GroupCat');
        pj = pairContrasts(m, pairs(k,:));
        if isnan(pLOO(k)) || pj > pLOO(k)
            pLOO(k) = pj; looMouse(k) = Tm.mouse(j);
        end
    end
end

fprintf('\n===== PRIMARY: per-mouse means, fitlm (df = %d), two-sided =====\n', mdl.DFE);
fprintf('%-22s %7s %8s %8s | %8s %8s %8s %8s %8s\n', 'Contrast', 'Delta', 'p', 'p Holm', ...
    'pair t', 'Welch', 'ranksum', 'median', 'LOO max');
for k = 1:nP
    fprintf('%-22s %+7.2f %8.4f %8.4f | %8.4f %8.4f %8.4f %8.4f %8.4f (drop %s)\n', ...
        [pairs{k,1} ' vs ' pairs{k,2}], dPri(k), pPri(k), pHolm(k), ...
        pT(k), pW(k), pRS(k), pMed(k), pLOO(k), char(looMouse(k)));
end

%% 4. For comparison: the previous ML fit on the same per-mouse table
lmeML = fitlme(Tm, 'meanDrift ~ GroupCat');
cn = lmeML.CoefficientNames; E = fixedEffects(lmeML);
cr = @(g) double(strcmp(cn, ['GroupCat_' g]));
fprintf('\n----- Previous fit (fitlme, ML, Satterthwaite) on the same table -----\n');
for k = 1:nP
    h = cr(pairs{k,2}) - cr(pairs{k,1});
    [p,~,~,df2] = coefTest(lmeML, h, 0, 'DFMethod', 'Satterthwaite');
    fprintf('  %-22s Delta = %+6.2f | df = %.1f | p = %.4f\n', ...
        [pairs{k,1} ' vs ' pairs{k,2}], h*E, df2, p);
end

%% 5. Field-level model (REML), for reference only
lmeF = fitlme(T, 'shiftData ~ GroupCat + (1|mouseID)', 'FitMethod', 'REML');
cnF = lmeF.CoefficientNames; EF = fixedEffects(lmeF);
crF = @(g) double(strcmp(cnF, ['GroupCat_' g]));
fprintf('\n----- Field-level LME (REML, Satterthwaite), reference only -----\n');
for k = 1:nP
    h = crF(pairs{k,2}) - crF(pairs{k,1});
    [p,~,~,df2] = coefTest(lmeF, h, 0, 'DFMethod', 'Satterthwaite');
    fprintf('  %-22s Delta = %+6.2f | df = %.0f | p = %.5f\n', ...
        [pairs{k,1} ' vs ' pairs{k,2}], h*EF, df2, p);
end

%% ===================== LOCAL FUNCTIONS =====================
function [p, d] = pairContrasts(mdl, pairs)
% Pairwise group contrasts from a fitlm model with reference level Control.
cn = mdl.CoefficientNames; b = mdl.Coefficients.Estimate;
cr = @(g) double(strcmp(cn, ['GroupCat_' g]));
p = nan(size(pairs,1),1); d = p;
for k = 1:size(pairs,1)
    h = cr(pairs{k,2}) - cr(pairs{k,1});
    p(k) = coefTest(mdl, h);
    d(k) = h * b;
end
end

function pa = holm(p)
% Holm step-down adjustment.
[ps, ix] = sort(p(:)); m = numel(ps);
adj = min(1, cummax((m:-1:1)' .* ps));
pa = nan(size(p)); pa(ix) = adj;
end