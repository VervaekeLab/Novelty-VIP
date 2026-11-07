%% ========================================================================
%% FIGURE 3b: timeline + single Novel-1 window contrast per group
%% ========================================================================
% STATS: collapse the four Novel-1 bins (x33-x48) to ONE value per mouse, then
% compare groups across mice (Corr ~ group). Two planned contrasts total
% (ArchT vs Control, ChrimsonR vs Control) -> no multiple-comparison question.
% PLOT: full timeline (mean +/- SEM across mice) with a Novel-1 window bracket
% carrying the two contrast p-values (no per-bin asterisks).
% Mouse IDs come from a REAL identifier column in tbl.
% ========================================================================
clearvars -except tbl; clc;

%% --- 1. DATA + REAL MOUSE ID --------------------------------------------
T = tbl; T.group = lower(string(T.group));
vn = T.Properties.VariableNames;
mouseCol = vn(contains(lower(vn),'mouse') | contains(lower(vn),'animal'));
if isempty(mouseCol), error('No mouse identifier column. Available: %s', strjoin(vn,', ')); end
T.MouseID = double(categorical(string(T.(mouseCol{1}))));
mice = unique(T.MouseID);
mGrp = strings(numel(mice),1);
for i = 1:numel(mice), mGrp(i) = T.group(find(T.MouseID==mice(i),1)); end
n1_bins = {'x33','x38','x43','x48'};

%% --- 2. NOVEL-1 WINDOW CONTRAST (one value per mouse, across mice) -------
pm = nan(numel(mice),1);
for i = 1:numel(mice)
    pm(i) = mean(mean(T{T.MouseID==mice(i), n1_bins},'omitnan'),'omitnan');
end
Tm = table(pm, categorical(mGrp,{'ctrl','arch','chri'},{'Control','ArchT','ChrimsonR'}), ...
           'VariableNames',{'Corr','GroupCat'});
Tm.GroupCat = reordercats(Tm.GroupCat,{'Control','ArchT','ChrimsonR'});
lme = fitlme(Tm,'Corr ~ GroupCat'); cn = lme.CoefficientNames; E = fixedEffects(lme);
DFM = {'DFMethod','Satterthwaite'};
P = struct();
fprintf('<strong>===== FIGURE 3B: Novel-1 window contrast (across mice) =====</strong>\n');
for op = {'ArchT','ChrimsonR'}
    h = double(strcmp(cn,['GroupCat_' op{1}]));
    [p,~,~,df2] = coefTest(lme,h,0,DFM{:}); est = h*E;
    P.(op{1}) = p;
    fprintf('  %-9s vs Control : est = %+6.3f | df = %.1f | p = %.4f\n', op{1}, est, df2, p);
end

%% --- 3. TIMELINE DATA (mean +/- SEM ACROSS MICE) ------------------------
allVars = T.Properties.VariableNames;
binCols = allVars(~cellfun(@isempty, regexp(allVars,'^x\d+$')));
x_coords = cellfun(@(c) str2double(c(2:end)), binCols);
[x_coords, sortIdx] = sort(x_coords); binCols = binCols(sortIdx);
names = {'Control','ArchT','ChrimsonR'}; keyw = {'ctrl','arch','chri'};
mu = zeros(3,numel(binCols)); se = zeros(3,numel(binCols));
for g = 1:3
    gm = mice(mGrp==keyw{g}); M = nan(numel(gm),numel(binCols));
    for i = 1:numel(gm), M(i,:) = mean(T{T.MouseID==gm(i), binCols},1,'omitnan'); end
    mu(g,:) = mean(M,1,'omitnan');
    se(g,:) = std(M,0,1,'omitnan') ./ sqrt(sum(~isnan(M),1));
end

%% --- 4. RENDER ----------------------------------------------------------
c_ctrl=[0 0 0]; c_arch=[0.18 0.45 0.71]; c_chri=[0.85 0.33 0.10];
colors = {c_ctrl,c_arch,c_chri};
fig = figure('Units','centimeters','Position',[5 5 12 13],'Color','w');
ax = axes('Parent',fig); hold(ax,'on');

% Novel-1 window shading (the tested window, x33-x48)
patch(ax,[31 50 50 31],[0 0 1.1 1.1],[0.90 0.90 0.90],'EdgeColor','none','HandleVisibility','off');

for g = 1:3
    errorbar(ax,x_coords,mu(g,:),se(g,:),'Color',colors{g},'LineWidth',2, ...
        'Marker','o','MarkerSize',7,'MarkerFaceColor','none','MarkerEdgeColor',colors{g}, ...
        'CapSize',5,'DisplayName',names{g});
end

% Novel-1 window bracket spanning the tested window (x33-x48), with a blue
% (ArchT) and a red (ChrimsonR) asterisk stacked vertically at the centre.
yb = 0.85;
plot(ax,[33 48],[yb yb],'-k','LineWidth',1.4,'HandleVisibility','off');
plot(ax,[33 33],[yb-0.02 yb],'-k','LineWidth',1.4,'HandleVisibility','off');
plot(ax,[48 48],[yb-0.02 yb],'-k','LineWidth',1.4,'HandleVisibility','off');
text(ax,40.5,yb+0.075,'*','Color',c_arch,'FontSize',24,'FontWeight','bold', ...
     'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','middle');
text(ax,40.5,yb+0.025,'*','Color',c_chri,'FontSize',24,'FontWeight','bold', ...
     'FontName','Arial','HorizontalAlignment','center','VerticalAlignment','middle');

title(ax,'Correlated to Novel 2','FontSize',16,'FontWeight','normal','FontName','Arial');
ylabel(ax,'Correlation coefficient','FontSize',16,'FontName','Arial');
xlabel(ax,'Trials (blocks of 5 trials)','FontSize',16,'FontName','Arial');
text(ax,15,1.05,'Familiar','FontSize',15,'FontName','Arial','HorizontalAlignment','center');
text(ax,40,1.05,'Novel 1','FontSize',15,'FontName','Arial','Color',[0.85 0.15 0.15],'HorizontalAlignment','center');
text(ax,60,1.05,'Novel 2','FontSize',15,'FontName','Arial','HorizontalAlignment','center');
legend(ax,names,'Location','southeast','Box','off','FontSize',12,'FontName','Arial');
set(ax,'TickDir','out','LineWidth',1.5,'FontSize',14,'FontName','Arial','Box','off', ...
       'XLim',[0 70],'YLim',[0 1.12],'YTick',0:0.2:1.0,'Layer','top');
hold(ax,'off');
