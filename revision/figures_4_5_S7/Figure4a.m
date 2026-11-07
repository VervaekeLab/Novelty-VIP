%% ========================================================================
%% FIGURE 3a: place-field emergence across the novel environment
%% ========================================================================
% Main panel : cumulative emergence curves, mean +/- SEM across mice.
% Novel 1    : level bracket, ArchT vs Control (delayed onset).
% Novel 2    : inset bar plot of the per-mouse emergence rate (slope), with the
%              ArchT-vs-Control slope contrast.
% Stats      : per-block linear mixed-effects models, evaluated across mice
%              (indTrCummData ~ GroupCat*trialC + (trialC|mouseID); Satterthwaite df).
% ========================================================================
clearvars -except lmeStats; clc;

%% --- 1. DATA PREP ---
T1=lmeStats.tables(1).tbl; T2=lmeStats.tables(2).tbl; T3=lmeStats.tables(3).tbl;
T=[T1;T2;T3];
rg=lower(string(T.group)); G=strings(height(T),1);
G(contains(rg,'arch'))="ArchT";
G(contains(rg,'chri'))="ChrimsonR";
G((contains(rg,'ctrl')|contains(rg,'control'))&G=="")="Control";
T.GroupCat=reordercats(categorical(G),{'Control','ArchT','ChrimsonR'});
T=T(~isundefined(T.GroupCat),:);
T.mouseID=categorical(T.mouseID);

c_ctrl=[0.20 0.20 0.20]; c_arch=[0.00 0.45 0.74]; c_chri=[0.85 0.33 0.10];
colors={c_ctrl,c_arch,c_chri}; names={'Control','ArchT','ChrimsonR'};

%% --- 2. CONTINUOUS CURVES (mouse-level: average within mouse, then across mice) ---
trials=1:60; cumm_mean=nan(60,3); cumm_sem=nan(60,3);
for t=trials
    for g=1:3
        sub = T(T.trialInd==t & T.GroupCat==names{g}, :);
        um  = unique(sub.mouseID);
        mv  = arrayfun(@(mm) mean(sub.indTrCummData(sub.mouseID==mm),'omitnan'), um);
        cumm_mean(t,g)=mean(mv,'omitnan');
        cumm_sem(t,g) =std(mv,'omitnan')/sqrt(numel(mv));
    end
end

%% --- 3. DATA-DRIVEN ANNOTATION STATS (per-context LME, Satterthwaite) ---
pN1_level = ctxContrast(T,[ 1 20],'level');   % ArchT vs control LEVEL in Novel 1
pN2_slope = ctxContrast(T,[21 40],'slope');   % ArchT vs control SLOPE in Novel 2
fprintf('N1: ArchT vs control level p = %.4f\n', pN1_level);
fprintf('N2: ArchT vs control slope p = %.4f\n', pN2_slope);

%% --- 3b. ACROSS-CONTEXT MODEL: group x context x trial interaction ---
% Single model across Novel 1-3 testing whether the trial-by-trial emergence
% rate depends on group and context (the omnibus interaction localized below).
A = T(T.trialInd<=60, :);                                  % Novel 1-3 (trials 1-60)
A.context  = categorical(ceil(double(A.trialInd)/20));     % 1,2,3 = Novel 1,2,3
A.trialC   = double(A.trialInd) - mean(double(A.trialInd),'omitnan');
A.GroupCat = removecats(A.GroupCat);
lmeFull = fitlme(A, 'indTrCummData ~ GroupCat*context*trialC + (trialC|mouseID:context)');
aov = anova(lmeFull,'DFMethod','Satterthwaite');
row = strcmp(aov.Term,'GroupCat:context:trialC');
fprintf('group x context x trial: F(%d,%.1f)=%.2f, p=%.4f\n', ...
        aov.DF1(row), aov.DF2(row), aov.FStat(row), aov.pValue(row));

%% --- 4. PER-MOUSE NOVEL-2 EMERGENCE SLOPES (the tested quantity; for the inset) ---
% Each mouse's own emergence rate across Novel 2 = slope of its cumulative
% curve over trials 21-40. These are the per-animal values underlying the
% LME slope contrast (pN2_slope), so the inset plots exactly what is tested.
n2 = 21:40; slopeByGroup = cell(1,3);
for g=1:3
    subG = T(T.GroupCat==names{g}, :);
    um   = unique(subG.mouseID);
    sl   = nan(numel(um),1);
    for i=1:numel(um)
        yv = arrayfun(@(t) mean(subG.indTrCummData(subG.mouseID==um(i) & subG.trialInd==t),'omitnan'), n2);
        ok = ~isnan(yv);
        if nnz(ok)>=2, pf = polyfit(n2(ok), yv(ok), 1); sl(i) = pf(1); end
    end
    slopeByGroup{g} = sl(~isnan(sl));
end

%% --- 5. FIGURE ---
fig=figure('Units','centimeters','Position',[5 5 16 16],'Color','w');
ax=axes('Parent',fig); hold(ax,'on');

% epoch shading (Novel 1, Novel 3)
patch(ax,[0 20.5 20.5 0],   [0 0 1.25 1.25],[0.94 0.94 0.94],'EdgeColor','none','HandleVisibility','off');
patch(ax,[40.5 60 60 40.5], [0 0 1.25 1.25],[0.94 0.94 0.94],'EdgeColor','none','HandleVisibility','off');

% curves + SEM bands
lh=gobjects(3,1);
for g=1:3
    mu=cumm_mean(:,g)'; se=cumm_sem(:,g)'; x=trials;
    fill(ax,[x fliplr(x)],[mu+se fliplr(mu-se)],colors{g}, ...
         'FaceAlpha',0.18,'EdgeColor','none','HandleVisibility','off');
    lh(g)=plot(ax,x,mu,'Color',colors{g},'LineWidth',2.5,'DisplayName',names{g});
end

% Novel-1 level bracket: ArchT vs Control
yb=1.06;
plot(ax,[2 19],[yb yb],'-','Color',c_arch,'LineWidth',1.6,'HandleVisibility','off');
plot(ax,[2 2],  [yb-0.018 yb],'-','Color',c_arch,'LineWidth',1.6,'HandleVisibility','off');
plot(ax,[19 19],[yb-0.018 yb],'-','Color',c_arch,'LineWidth',1.6,'HandleVisibility','off');
text(ax,14,yb+0.012,sprintf('ArchT level, p = %.3f',pN1_level),'Color',c_arch, ...
     'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',14,'FontName','Arial');

% epoch headers
text(ax,10,1.21,'Novel 1','HorizontalAlignment','center','FontSize',14,'FontWeight','bold','FontName','Arial','Color',[0.85 0.15 0.15]);
text(ax,30,1.21,'Novel 2','HorizontalAlignment','center','FontSize',14,'FontWeight','bold','FontName','Arial','Color','k');
text(ax,50,1.21,'Novel 3','HorizontalAlignment','center','FontSize',14,'FontWeight','bold','FontName','Arial','Color',[0.85 0.15 0.15]);

ylabel(ax,'Emergence of new place fields','FontSize',15,'FontName','Arial','FontWeight','bold');
xlabel(ax,'Trials','FontSize',15,'FontName','Arial','FontWeight','bold');
ylim(ax,[0 1.25]); xlim(ax,[0 60]); yticks(ax,0:0.2:1.0);
% tight colour-coded text legend (no line samples), lower area of Novel 1/2
text(ax,13,0.27,'Control',  'Color',c_ctrl,'FontWeight','bold','FontSize',13,'FontName','Arial','HorizontalAlignment','left');
text(ax,13,0.18,'ArchT',    'Color',c_arch,'FontWeight','bold','FontSize',13,'FontName','Arial','HorizontalAlignment','left');
text(ax,13,0.09,'ChrimsonR','Color',c_chri,'FontWeight','bold','FontSize',13,'FontName','Arial','HorizontalAlignment','left');
set(ax,'TickDir','out','LineWidth',1.5,'FontSize',13,'Layer','top','FontName','Arial','Box','off');
hold(ax,'off');

%% --- Novel-2 emergence-rate INSET (bars = mean +/- SEM, points = mice) ---
axI = axes('Parent',fig,'Position',[0.59 0.180 0.275 0.360]); hold(axI,'on');
gl = {'Ctrl','ArchT','Chr'}; allSl = [];
for g=1:3
    sl = slopeByGroup{g}; allSl = [allSl; sl(:)];
    mu = mean(sl); se = std(sl)/sqrt(numel(sl));
    bar(axI,g,mu,0.62,'FaceColor',colors{g},'FaceAlpha',0.30, ...
        'EdgeColor',colors{g},'LineWidth',1.1);
    plot(axI,[g g],[mu mu+se],'-','Color',colors{g},'LineWidth',1.2);          % SEM whisker
    plot(axI,[g-0.10 g+0.10],[mu+se mu+se],'-','Color',colors{g},'LineWidth',1.2);
    xj = g + (rand(numel(sl),1)-0.5)*0.20;
    scatter(axI,xj,sl,20,colors{g},'filled','MarkerEdgeColor','w','LineWidth',0.4);
end
yhi = max(allSl); pad = 0.10*yhi;
plot(axI,[1 2],[yhi+pad yhi+pad],'-k','LineWidth',1);   % ArchT vs control
text(axI,1.5,yhi+1.8*pad,sprintf('p = %.3f',pN2_slope), ...
     'HorizontalAlignment','center','FontSize',12,'FontName','Arial');
xlim(axI,[0.4 3.6]); ylim(axI,[0 yhi+3.4*pad]);
xticks(axI,1:3); xticklabels(axI,gl);
ylabel(axI,{'Novel 2','place cells/trial'},'FontSize',12,'FontName','Arial');
set(axI,'TickDir','out','LineWidth',1,'FontSize',11,'FontName','Arial','Box','off');

%% ===================== LOCAL FUNCTIONS =====================
function p = ctxContrast(T, trange, kind)
% ArchT-vs-control level or slope contrast within one context, Satterthwaite.
    S = T(T.trialInd>=trange(1) & T.trialInd<=trange(2), :);
    S.trialC   = double(S.trialInd) - mean(double(S.trialInd),'omitnan');  % center within context
    S.GroupCat = removecats(S.GroupCat);
    lme = fitlme(S, 'indTrCummData ~ GroupCat*trialC + (trialC|mouseID)');
    cn  = lme.CoefficientNames;
    if strcmp(kind,'level'), tok={'GroupCat_ArchT'}; else, tok={'GroupCat_ArchT','trialC'}; end
    p = coefTest(lme, coefRow(cn,tok), 0, 'DFMethod','Satterthwaite');
end

function v = coefRow(cn, tokens)
% Order-independent coefficient selector (handles a:b vs b:a naming).
    tokens=string(tokens); v=zeros(1,numel(cn));
    for k=1:numel(cn)
        parts=string(strsplit(cn{k},':'));
        if numel(parts)==numel(tokens) && isempty(setdiff(parts,tokens)) ...
                                       && isempty(setdiff(tokens,parts))
            v(k)=1; return;
        end
    end
end
