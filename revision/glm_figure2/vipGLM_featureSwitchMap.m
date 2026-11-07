function vipGLM_featureSwitchMap()
    % default folders: cohort RESULTS in the MATLAB userpath (.../Documents/MATLAB),
    % session DATA in Desktop/VIP-GCaMP-sData.
    dfRes=userpath; if isempty(dfRes)||~exist(dfRes,'dir'), dfRes=pwd; end
    dfSess=getenv('HOME'); if isempty(dfSess), dfSess=getenv('USERPROFILE'); end
    dfSess=fullfile(dfSess,'Desktop','VIP-GCaMP-sData'); if ~exist(dfSess,'dir'), dfSess=dfRes; end
% VIPGLM_FEATURESWITCHMAP  Two stacked heatmaps (features x cells) of unique
% dR2, FAMILIAR on top and NOVEL directly below, sharing the same column
% order, so you can read down each column to see whether a cell's encoding
% profile changes from familiar to novel.
%
% INPUT  the COHORT RESULTS file (ONE .mat). OUTPUT two stacked hetmaps.
%
% Columns = cells, sorted by FAMILIAR primary feature (grouped).
% Rows (7) = Vel(+), Vel(-), Accel, Licking, Reward, Position, Pupil.
%   Velocity dR2 is placed in the Vel(+) or Vel(-) row per the cell's velocity
%   sign (saved velSign; falls back to the signed primary label if absent).
%
% Input (prompted): saved cohort results .mat with res.cellTable.
% No re-run / no session files needed.

        uiwait(msgbox({'Select the COHORT RESULTS file.','','ONE file: the cohortGLM_*.mat saved by vipGLM_cohort (has res.cellTable).','','NOT a session/m8xxx file.'},'Cohort results (ONE file)','help'));
    [rf,rp]=uigetfile('*.mat','COHORT RESULTS file (ONE .mat, cohortGLM_*.mat)',dfRes);
    if isequal(rf,0), return; end
    R=load(fullfile(rp,rf));
    if ~isfield(R,'res')||~isfield(R.res,'cellTable'), error('No res.cellTable.'); end
    T=R.res.cellTable;

    featCol={'dR2_Velocity','dR2_Acceleration','dR2_Licking','dR2_Reward','dR2_Position','dR2_Pupil'};
    for k=1:numel(featCol)
        if ~any(strcmpi(T.Properties.VariableNames,featCol{k}))
            error('cellTable missing column %s', featCol{k});
        end
    end
    rowNames={'Vel (+)','Vel (-)','Accel','Licking','Reward','Position','Pupil'};
    nRow=7;

    Tf=T(T.context=="Familiar" & T.ok,:);
    Tn=T(T.context=="Novel"    & T.ok,:);
    % match cells present in BOTH contexts by (session,cellID)
    keyF=string(Tf.session)+"|"+string(Tf.cellID);
    keyN=string(Tn.session)+"|"+string(Tn.cellID);
    [common,iaF,~]=intersect(keyF,keyN,'stable');
    if isempty(common), error('No cells present in both contexts.'); end
    Tf=Tf(iaF,:);
    % reorder Tn to match Tf order
    [~,locN]=ismember(keyF(iaF),keyN); Tn=Tn(locN,:);
    nCells=height(Tf);

    Mf=buildMat(Tf,featCol,nRow);   % 7 x nCells familiar
    Mn=buildMat(Tn,featCol,nRow);   % 7 x nCells novel (same cell order as Tf for now)

    % sort columns by familiar primary feature (row index of max familiar dR2),
    % then by descending that max within each group
    [~,primRow]=max(Mf,[],1);
    peak=max(Mf,[],1);
    [~,order]=sortrows([primRow(:), -peak(:)]);
    Mf=Mf(:,order); Mn=Mn(:,order); primRowSorted=primRow(order);

    % shared color scale across both matrices
    cmax=prctile([Mf(:);Mn(:)],99); if cmax<=0, cmax=max([Mf(:);Mn(:);eps]); end

    figure('Color','w','Position',[60 80 1180 620]);
    ax1=subplot(2,1,1);
    imagesc(Mf,[0 cmax]); colormap(ax1,parula);
    set(ax1,'YTick',1:nRow,'YTickLabel',rowNames,'TickDir','out');
    title('Familiar: unique \DeltaR^2 (rows=features, cols=cells sorted by familiar primary)');
    xlabel('VIP cell'); cb1=colorbar; cb1.Label.String='\DeltaR^2';
    drawGroupDividers(primRowSorted, nRow);

    ax2=subplot(2,1,2);
    imagesc(Mn,[0 cmax]); colormap(ax2,parula);
    set(ax2,'YTick',1:nRow,'YTickLabel',rowNames,'TickDir','out');
    title('Novel: unique \DeltaR^2 (same cell order as above)');
    xlabel('VIP cell'); cb2=colorbar; cb2.Label.String='\DeltaR^2';
    drawGroupDividers(primRowSorted, nRow);

    sgtitle(sprintf('Per-cell feature encoding, familiar vs novel (n=%d cells in both contexts)',nCells), ...
        'FontWeight','bold');
end

function M = buildMat(Tx, featCol, nRow)
    n=height(Tx); M=zeros(nRow,n);
    getc=@(name) double(Tx.(name));
    vel=getc(featCol{1}); acc=getc(featCol{2}); lick=getc(featCol{3});
    rew=getc(featCol{4}); pos=getc(featCol{5}); pup=getc(featCol{6});
    % velocity sign -> route to row 1 (Vel+) or row 2 (Vel-)
    vs=velSignVec(Tx);
    for c=1:n
        if vs(c)>=0, M(1,c)=vel(c); else, M(2,c)=vel(c); end
        M(3,c)=acc(c); M(4,c)=lick(c); M(5,c)=rew(c); M(6,c)=pos(c); M(7,c)=pup(c);
    end
    M(~isfinite(M))=0; M(M<0)=0;   % clamp tiny negative dR2 to 0 for display
end

function vs = velSignVec(Tx)
    n=height(Tx); vs=ones(n,1);
    if any(strcmpi(Tx.Properties.VariableNames,'velSign'))
        v=double(Tx.velSign); v(isnan(v))=1; vs=sign(v); vs(vs==0)=1; return;
    end
    % fallback: parse signed primary label if present
    if any(strcmpi(Tx.Properties.VariableNames,'primaryFeatureSigned'))
        s=string(Tx.primaryFeatureSigned);
        vs(contains(s,"(-)"))=-1;
    end
end

function drawGroupDividers(primRowSorted, nRow)
    hold on;
    b=find(diff(primRowSorted(:))~=0);   % column boundaries between groups
    for k=1:numel(b)
        x=b(k)+0.5; plot([x x],[0.5 nRow+0.5],'w-','LineWidth',1.5);
    end
    hold off;
end
