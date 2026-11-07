function vipGLM_primarySecondary(matFile)
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_PRIMARYSECONDARY  Post-hoc primary/secondary feature analysis.
%   Operates on the saved cohort results (res.cellTable) -- NO re-running.
%
%   INPUT  the COHORT RESULTS file: pass its path as matFile, or call with no
%          args to get a file prompt (ONE .mat). OUTPUT printed analysis.
%
%   Answers two reciprocal questions, within each context:
%     (1) For cells whose PRIMARY feature is Pupil, what is their SECONDARY?
%         (Is it Velocity? -> would reflect pupil/velocity collinearity.)
%     (2) For cells whose PRIMARY feature is Velocity, is Pupil the SECONDARY?
%
%   "Secondary" = the significant feature with the 2nd-largest unique dR2
%   (primary = largest). Cells with only one significant feature have no
%   secondary and are reported separately.
%
%   Usage:
%     vipGLM_primarySecondary                 % prompts for the .mat file
%     vipGLM_primarySecondary('cohortGLM_10sessions_...mat')

    if nargin<1 || isempty(matFile)
            uiwait(msgbox({'Select the COHORT RESULTS file.','','ONE file: the cohortGLM_*.mat saved by vipGLM_cohort (has res.cellTable).','','NOT a session/m8xxx file.'},'Cohort results (ONE file)','help'));
        [f,p]=uigetfile('*.mat','COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes);
        if isequal(f,0), return; end
        matFile=fullfile(p,f);
    end
    S=load(matFile);
    if ~isfield(S,'res') || ~isfield(S.res,'cellTable')
        error('File does not contain res.cellTable. Run vipGLM_cohort with Save=Yes first.');
    end
    T=S.res.cellTable;
    feats={'Velocity','Acceleration','Licking','Reward','Position','Pupil'};

    fprintf('\n==================================================================\n');
    fprintf('Primary / secondary feature analysis\n  %s\n', matFile);
    fprintf('==================================================================\n');

    contexts=unique(T.context,'stable');
    for ci=1:numel(contexts)
        ctx=contexts(ci);
        Tc=T(T.context==ctx & T.ok, :);
        fprintf('\n----- Context: %s  (%d fitted cells) -----\n', ctx, height(Tc));

        % compute primary & secondary feature per cell from the dR2/sig columns
        [primF, secF] = deal(strings(height(Tc),1));
        for i=1:height(Tc)
            [primF(i), secF(i)] = primSec(Tc(i,:), feats);
        end

        % (1) Pupil-primary -> secondary distribution
        reportDirection('Pupil', primF, secF, feats);
        % (2) Velocity-primary -> secondary distribution
        reportDirection('Velocity', primF, secF, feats);

        % focused reciprocal summary
        nPupVel = sum(primF=="Pupil"   & secF=="Velocity");
        nPup    = sum(primF=="Pupil");
        nVelPup = sum(primF=="Velocity"& secF=="Pupil");
        nVel    = sum(primF=="Velocity");
        fprintf('\n  >> Reciprocal summary (%s):\n', ctx);
        fprintf('     Pupil-primary with Velocity secondary: %d / %d (%.0f%%)\n', ...
            nPupVel, nPup, 100*safeDiv(nPupVel,nPup));
        fprintf('     Velocity-primary with Pupil secondary: %d / %d (%.0f%%)\n', ...
            nVelPup, nVel, 100*safeDiv(nVelPup,nVel));
    end
    fprintf('\nDone.\n');
end

% ---- primary & secondary feature for one cell row ----
% Uses unique dR2 ranked among SIGNIFICANT features. Returns "none" where
% there is no primary / no secondary.
function [prim, sec] = primSec(row, feats)
    prim="none"; sec="none";
    nF=numel(feats); dr=nan(1,nF); sg=false(1,nF);
    for f=1:nF
        dr(f)=row.(['dR2_' feats{f}]);
        sg(f)=logical(row.(['sig_' feats{f}]));
    end
    valid=find(sg & ~isnan(dr));
    if isempty(valid), return; end
    [~,order]=sort(dr(valid),'descend');
    ranked=valid(order);
    prim=string(feats{ranked(1)});
    if numel(ranked)>=2, sec=string(feats{ranked(2)}); end
end

function reportDirection(primaryName, primF, secF, feats)
    sel = primF==primaryName;
    nTot=sum(sel);
    fprintf('\n  Cells with PRIMARY = %s: n=%d\n', primaryName, nTot);
    if nTot==0, return; end
    secOfThese=secF(sel);
    nNone=sum(secOfThese=="none");
    fprintf('    secondary feature breakdown:\n');
    for f=1:numel(feats)
        c=sum(secOfThese==feats{f});
        if c>0
            fprintf('      %-14s %3d  (%.0f%% of %s-primary)\n', feats{f}, c, 100*c/nTot, primaryName);
        end
    end
    if nNone>0
        fprintf('      %-14s %3d  (%.0f%%) -- only one significant feature\n', '(none)', nNone, 100*nNone/nTot);
    end
end

function r = safeDiv(a,b), if b==0, r=0; else, r=a/b; end, end
