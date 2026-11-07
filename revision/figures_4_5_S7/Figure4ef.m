%% ========================================================================
%% FIGURE 4e, f: per-block speed (4e) and lick (4f) group analysis
%% ========================================================================
% Question : in each trial block, does ArchT silencing or ChrimsonR excitation
%            change reward-predictive speed (4e) / lick (4f) vs Control?
% Model    : fit separately per block and per behaviour, evaluated across mice
%            y ~ GroupCat + trialC + (1|mouseID); Satterthwaite df.
%              trialC      : within-block trial index, mean-centred (covariate)
%              (1|mouseID) : animal random intercept
% Contrasts: planned, pre-specified; exact p-values; no family-wise correction.
%              ArchT vs Control     : one-sided (pre-specified decrease)
%              ChrimsonR vs Control : two-sided
% Data     : lmeStats.tables(1:5).tbl = Fam1, Novel1, Novel2, Novel3, Fam2.
% ========================================================================
clearvars -except lmeStats; clc;

blockNames = {'Fam 1','Novel 1','Novel 2','Novel 3','Fam 2'};
behaviors  = {'speedData','lickData'};
behLabels  = {'Speed','Lick'};
archDir    = -1;   % ArchT pre-specified direction (decrease)

% stored as [block x behaviour x (1 = ArchT one-sided, 2 = ChrimsonR two-sided)]
con_est = NaN(5,2,2); con_p = NaN(5,2,2);

fprintf('<strong>============== FIGURE 4: PER-BLOCK GROUP CONTRASTS ==============</strong>\n');
fprintf('y ~ GroupCat + trialC + (1|mouseID); Satterthwaite df; exact p; no correction.\n');
fprintf('ArchT vs Control: one-sided (pre-specified decrease).  ChrimsonR vs Control: two-sided.\n');

for k = 1:numel(blockNames)
    T = lmeStats.tables(k).tbl;
    rg = lower(string(T.group)); G = strings(height(T),1);
    G(contains(rg,'ctrl')|contains(rg,'control')) = "Control";
    G(contains(rg,'arch')) = "ArchT";
    G(contains(rg,'chri')) = "ChrimsonR";
    T.GroupCat = reordercats(categorical(G), {'Control','ArchT','ChrimsonR'});
    T.mouseID  = categorical(string(T.mouseID));
    T.trialC   = double(T.trialInd) - mean(double(T.trialInd),'omitnan');

    fprintf('\n----- BLOCK %d: %s -----\n', k, blockNames{k});
    for v = 1:2
        y   = behaviors{v};
        lme = fitlme(T, sprintf('%s ~ GroupCat + trialC + (1|mouseID)', y));
        cn  = lme.CoefficientNames; E = fixedEffects(lme); DFM = {'DFMethod','Satterthwaite'};

        % ArchT vs Control: one-sided (pre-specified decrease)
        hA = coefRow(cn,{'GroupCat_ArchT'}); eA = hA*E;
        p2 = coefTest(lme, hA, 0, DFM{:});                  % two-sided
        if sign(eA)==archDir, pA = p2/2; else, pA = 1 - p2/2; end

        % ChrimsonR vs Control: two-sided
        hC = coefRow(cn,{'GroupCat_ChrimsonR'}); eC = hC*E;
        pC = coefTest(lme, hC, 0, DFM{:});

        con_est(k,v,1)=eA; con_p(k,v,1)=pA;
        con_est(k,v,2)=eC; con_p(k,v,2)=pC;

        fprintf('  %-6s | ArchT %+7.3f  p=%.4f (1-sided) | ChrimsonR %+7.3f  p=%.4f (2-sided)\n', ...
            behLabels{v}, eA, pA, eC, pC);
    end
end

fprintf('\n================================================================\n');
fprintf('Stored: con_est, con_p  [block x behaviour x (1=ArchT 1-sided, 2=ChrimsonR 2-sided)].\n');

%% ===================== LOCAL FUNCTION =====================
function v = coefRow(cn, tokens)
% Order-independent coefficient selector (robust to a:b vs b:a naming).
    tokens = string(tokens); v = zeros(1, numel(cn));
    for k = 1:numel(cn)
        parts = string(strsplit(cn{k}, ':'));
        if numel(parts)==numel(tokens) && isempty(setdiff(parts,tokens)) ...
                                       && isempty(setdiff(tokens,parts))
            v(k) = 1; return;
        end
    end
end
