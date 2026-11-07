%% FIGURE 3e/3f - PLACE-FIELD DRIFT GROUP ANALYSIS
% Data: tbl = native table in workspace containing:
%   shiftData (drift, cm) | group (categorical) | mouseID (categorical)
%
% DESIGN (animal-level, mixed-effects framework; pairwise contrasts only):
%   Group is a BETWEEN-mouse factor, so the comparison has at most ~n_mice
%   degrees of freedom. To keep the test at the animal level (and avoid the
%   field-level df inflation that a field-level LME produces for a between-mouse
%   factor), we fit the linear mixed-effects model on the per-mouse summaries:
%       meanDrift ~ GroupCat
%   This is the animal-level instance of the same mixed-model framework used
%   throughout (one summary per animal -> contrasts evaluated across mice).
%   Field-level LME on all place fields is reported separately as corroboration.

clearvars -except tbl; clc;

% --- 1. SETUP & DATA CLEANING ---
T = tbl;
gs = lower(string(T.group)); G = strings(height(T),1);
G(contains(gs,'ctrl')) = "Control";
G(contains(gs,'arch')) = "ArchT";
G(contains(gs,'chri')) = "ChrimsonR";
T.GroupCat = reordercats(categorical(G), {'Control','ArchT','ChrimsonR'});
T.mouseID  = categorical(T.mouseID);
T.logShift = log1p(T.shiftData);

% --- 2. PER-MOUSE SUMMARIES (the animal-level data; n must be 15) ---
mice = categories(removecats(T.mouseID));
pmDrift = nan(numel(mice),1); pmGrp = strings(numel(mice),1);
for i = 1:numel(mice)
    s = T(T.mouseID == mice{i}, :);
    pmDrift(i) = mean(s.shiftData,'omitnan');     % one value per mouse
    pmGrp(i)   = string(s.GroupCat(1));
end
fprintf('Per-mouse summaries: n = %d mice (Ctrl/Arch/Chr = %d/%d/%d)\n', numel(mice), ...
    sum(pmGrp=="Control"), sum(pmGrp=="ArchT"), sum(pmGrp=="ChrimsonR"));

Tm = table(pmDrift, reordercats(categorical(pmGrp),{'Control','ArchT','ChrimsonR'}), ...
           'VariableNames', {'meanDrift','GroupCat'});
pairs = {'Control','ArchT'; 'Control','ChrimsonR'; 'ArchT','ChrimsonR'};

%% ===== PRIMARY: ANIMAL-LEVEL MIXED-EFFECTS (per-mouse means, n=15) =====
% With one observation per mouse this is a linear model in the LME framework;
% the between-mouse contrast df is ~ n_mice - n_groups (NOT field-level).
fprintf('\n<strong>===== PRIMARY: ANIMAL-LEVEL LME (per-mouse means) =====</strong>\n');
lme = fitlme(Tm, 'meanDrift ~ GroupCat');
cn = lme.CoefficientNames; E = fixedEffects(lme); DFM = {'DFMethod','Satterthwaite'};
cr = @(g) double(strcmp(cn, ['GroupCat_' g]));
for k = 1:3
    h = cr(pairs{k,2}) - cr(pairs{k,1});
    [p,F,df1,df2] = coefTest(lme, h, 0, DFM{:});
    est = h*E;
    fprintf('  %-9s vs %-9s : Delta = %+6.2f cm | df = %.1f | p = %.4f %s\n', ...
        pairs{k,1}, pairs{k,2}, est, df2, p, star(p));
end

%% ===== CORROBORATION: FIELD-LEVEL LME (all place fields) =====
% Reported for completeness; note its contrast df is field-level (not animal-level).
fprintf('\n<strong>===== CORROBORATION: FIELD-LEVEL LME (Satterthwaite) =====</strong>\n');
lmeF = fitlme(T, 'shiftData ~ GroupCat + (1|mouseID)');
cnF = lmeF.CoefficientNames; EF = fixedEffects(lmeF);
crF = @(g) double(strcmp(cnF, ['GroupCat_' g]));
for k = 1:3
    h = crF(pairs{k,2}) - crF(pairs{k,1});
    [p,F,df1,df2] = coefTest(lmeF, h, 0, DFM{:});
    est = h*EF;
    fprintf('  %-9s vs %-9s : Delta = %+6.2f cm | df = %.0f | p = %.5f %s\n', ...
        pairs{k,1}, pairs{k,2}, est, df2, p, star(p));
end
fprintf('=================================================================\n');

%% ===================== LOCAL FUNCTION =====================
function s = star(p)
    if isnan(p), s=''; else, s = repmat('*',1,(p<0.05)+(p<0.01)+(p<0.001)); end
end
