function basic(do_correlation, do_scatter_plots, do_error_bars, do_anova)

    if nargin < 1
        do_correlation = false;
    end
    if nargin < 2
        do_scatter_plots = false;
    end
    if nargin < 3
        do_error_bars = true;
    end
    if nargin < 4
        do_anova = false;
    end

    T = importfile('data.csv');
    T.sp = categorical(T.sp);
    T.sun_shade = categorical(T.sun_shade);
    T.GroupHscon = categorical(T.GroupHscon);
    T.HsconRang = max(T.HsconRang) - T.HsconRang;
    
    %% SPECIFY OUTPUT FOLDERS
    
    barplot_dir = 'fig/bar/';
    scatter_dir = 'fig/scatter/';
    stat_dir = 'stat/';
    mkdir(stat_dir);
    mkdir(barplot_dir);
    mkdir(scatter_dir);
    global dark_green light_green vars names
    dark_green = [1 103 52]./255;
    light_green = [113 183 53]./255;
    
    %% SELECTE SUB TABLES
    
    [uniqueSp, uniqueIdx] = unique(T.sp);
    
    vars = {'Blattdicke_mm', 'Chlorophyll', 'Reissfestigkeit_N', ...
        'DW_FW', 'SLA', 'Stomatadichte', 'd15N14N', 'd13C12C', ...
        'N', 'C', 'PARsat', 'ETR_1500', 'Hzuwachs', 'init_slope'};
    names = {'Blattdicke ($mm$)', 'Chlorophyll (SPAD)', 'Reissfestigkeit ($N$)', ...
        'LDMC ($mg/g$)', 'SLA  ($mm^2/mg$)', 'Stomatadichte ($n/mm^2$)', ...
        '$d~15N/14N$', '$d~13C/12C$', 'N ($\%$)', 'C ($\%$)', 'PARsat', ...
        'ETR1500 ($\mu$mol $m^{-2}~s^{-1}$)', 'H\"ohenzuwachs ($cm/year$)', ...
        'Initial Slope'};
    varSel = ismember(T.Properties.VariableNames, vars);
    
    U = arrayfun(@(sp) T(T.sp == sp,:), uniqueSp, 'uniform', 0);
    TSun = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'sonne',:), ...
        uniqueSp, 'uniform', 0);
    TShade = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'schatten',:), ...
        uniqueSp, 'uniform', 0);
    
    [countSun,~] = cellfun(@size, TSun);
    [countShade,~] = cellfun(@size, TShade);
    
    sel = countSun >= 0 & countShade >= 0;
    selectedT = restrict(T, 0, 0);
    selectedSp = uniqueSp(sel);
    selectedSpGroup = T.HsconRang(uniqueIdx);
    selectedSpGroup = selectedSpGroup(sel);
    selectedSpHzuwachs = T.Hzuwachs(uniqueIdx);
    selectedSpHzuwachs = selectedSpHzuwachs(sel);
    selectedTSun = TSun(sel);
    selectedTShade = TShade(sel);
    
    meanTSun = cellfun(@(X) nanmean(table2array(X(:,vars)), 1), selectedTSun, 'uniform', 0);
    meanTShade = cellfun(@(X) nanmean(table2array(X(:,vars)), 1), selectedTShade, 'uniform', 0);
    meanTSun = vertcat(meanTSun{:});
    meanTShade = vertcat(meanTShade{:});
    
    selTSun = vertcat(selectedTSun{:});
    selTShade = vertcat(selectedTShade{:});
    
    %% PEARSON CORRELATION COEFFICIENTS
    
    if do_correlation
%         TSun = table2array(selTSun(:,vars));
        TSun = meanTSun;
        [X,Y] = meshgrid(1:numel(vars));
        [RSun,PSun] = arrayfun(@(x,y) corr(TSun(:,x),TSun(:,y), 'rows', 'complete'), X, Y);
        
%         TShade = table2array(selTShade(:,vars));
        TShade = meanTShade;
        [RShade,PShade] = arrayfun(@(x,y) corr(TShade(:,x),TShade(:,y), 'rows', 'complete'), X, Y);
        
        CorrSun = RSun;
        TF = triu(ones(size(PSun))) == 1;
        CorrSun(TF) = PSun(TF);
        CorrSun = addnames(array2table(CorrSun), names, vars);
        CorrShade = RShade;
        CorrShade(TF) = PShade(TF);
        CorrShade = addnames(array2table(CorrShade), names, vars);
        writetable(CorrSun, [stat_dir 'Corr.xlsx'], 'WriteRowNames',true, 'Sheet', 'Sun');
        writetable(CorrShade, [stat_dir 'Corr.xlsx'], 'WriteRowNames',true, 'Sheet', 'Shade');
    end
    
    %% SCATTER PLOT COMPARISONS
    
    if do_scatter_plots
%         TSunSp = selTSun(:,'sp');
%         TShadeSp = selTShade(:,'sp');
        
        pvalue = 0.05;
        saveAllSignificantCorr(TSun, [], TShade, [], CorrSun, CorrShade, pvalue, scatter_dir);
        close all;
    end
    
    %% SCATTER PLOT MATRIX
    
%     if 0
%         plotandsavematrix(selTSun(:,vars), names);
%         saveFigure([scatter_dir 'scattermatrix_sun'], [1000 1000], '-dpdf');
%         plotandsavematrix(selTShade(:,vars), names);
%         saveFigure([scatter_dir 'scattermatrix_shade'], [1000 1000], '-dpdf');
%         
%         close all
%     end
    
    %% PERFORM TTEST
    
    if do_error_bars
        compute_and_save_errorbars(T, barplot_dir);
%         nSun = cell2mat(cellfun(@(t) sum(~isnan(table2array(t(:,varSel))), 1), ...
%             selectedTSun, 'uniform', 0));
%         nShade = cell2mat(cellfun(@(t) sum(~isnan(table2array(t(:,varSel))), 1), ...
%             selectedTShade, 'uniform', 0));
% 
%         meanSun = cell2mat(cellfun(@(t) nanmean(table2array(t(:,varSel)), 1), ...
%             selectedTSun, 'uniform', 0));
%         meanShade = cell2mat(cellfun(@(t) nanmean(table2array(t(:,varSel)), 1), ...
%             selectedTShade, 'uniform', 0));
% 
%         stdSun = cell2mat(cellfun(@(t) nanstd(table2array(t(:,varSel)), 0, 1), ...
%             selectedTSun, 'uniform', 0));
%         stdShade = cell2mat(cellfun(@(t) nanstd(table2array(t(:,varSel)), 0, 1), ...
%             selectedTShade, 'uniform', 0));
%         
%         % sort by group
%         %sortIdx = 1:numel(selectedSpGroup);
%         [~,sortIdx] = sort(selectedSpGroup);
%         selectedSp = selectedSp(sortIdx);
%         meanSun = meanSun(sortIdx,:);
%         meanShade = meanShade(sortIdx,:);
%         stdSun = stdSun(sortIdx,:);
%         stdShade = stdShade(sortIdx,:);
%         nSun = nSun(sortIdx,:);
%         nShade = nShade(sortIdx,:);
% 
%         [h,p] = cellfun(@(u,h) ttest2(table2array(u(:,varSel)),table2array(h(:,varSel))), ...
%             selectedTSun(sortIdx), selectedTShade(sortIdx), 'uniform', 0);
%         p = cell2mat(p);
%         labels = cell(size(p));
%         labels(p < 0.05) = {'\textbf{*}'};
%         labels(p < 0.01) = {'\textbf{**}'};
%         labels(p < 0.001) = {'\textbf{***}'};
%         
%         for i = 1:size(meanSun,2)
%             x_labels = arrayfun(@(sp,l) [l{1} '\textsf{' char(sp) '}'], ...
%                 selectedSp, labels(:,i), 'uniform', 0);
%             errbarplot([meanSun(:,i) meanShade(:,i)], ...
%                 [stdSun(:,i) stdShade(:,i)], ...
%                 names{i}, x_labels);
%             saveFigure([barplot_dir vars{i}], [600 400], '-dsvg');
%         end
% 
%         idx = (1:size(meanSun,2)) * 2 - 1;
% 
%         mu = zeros(size(meanSun) .* [1 2]);
%         mu(:,idx) = meanSun;
%         mu(:,idx + 1) = meanShade;
%         mu = array2table(mu);
%         mu.Properties.VariableNames(idx) = cellfun(@(i) ['mu_sun_',i], vars, 'uniform', 0);
%         mu.Properties.VariableNames(idx + 1) = cellfun(@(i) ['mu_shade_',i], vars, 'uniform', 0);
% 
%         sd = zeros(size(stdSun) .* [1 2]);
%         sd(:,idx) = stdSun;
%         sd(:,idx + 1) = stdShade;
%         sd = array2table(sd);
%         sd.Properties.VariableNames(idx) = cellfun(@(i) ['sd_sun_',i], vars, 'uniform', 0);
%         sd.Properties.VariableNames(idx + 1) = cellfun(@(i) ['sd_shade_',i], vars, 'uniform', 0);
% 
%         num = zeros(size(stdSun) .* [1 2]);
%         num(:,idx) = nSun;
%         num(:,idx + 1) = nShade;
%         num = array2table(num);
%         num.Properties.VariableNames(idx) = cellfun(@(i) ['num_sun_',i], vars, 'uniform', 0);
%         num.Properties.VariableNames(idx + 1) = cellfun(@(i) ['num_shade_',i], vars, 'uniform', 0);
% 
%         p = array2table(p);
%         p.Properties.VariableNames = cellfun(@(i) ['p_',i], vars, 'uniform', 0);
% 
%         h = cell2mat(h);
%         h = array2table(h);
%         h.Properties.VariableNames = cellfun(@(i) ['h_',i], vars, 'uniform', 0);
% 
%         summary = [array2table(selectedSp), mu, sd, num, p, h];
%         writetable(summary, [stat_dir 'mean_std_ttest.xlsx'], 'WriteRowNames',true);
%         
%         close all
    end
    
    %% PERFORM ANOVA per trait
    
    if do_anova
        compute_and_save_anovas(T, stat_dir);
    end
    
    %% PERFORM ANCOVA
    
    if 0
        pvalue = 0.05;
        U = table2array(selectedT(:,vars));
%         G = categorical(arrayfun(...
%             @(a,b) [char(a) '-' char(b)], ...
%             selectedT.sp, selectedT.sun_shade, 'uniform', 0));
        G = selectedT.sun_shade;
        
        [X,Y] = meshgrid(1:numel(vars));
%         TF = ~(X == Y);
%         X = X(TF);
%         Y = Y(TF);
        [~,atab,~,~] = arrayfun(...
            @(x,y) aoctool(U(:,x),U(:,y),G,pvalue,vars{x},vars{y},'Sun/Shade','off'), ...
            X, Y, 'uniform', 0);
        F = cellfun(@(a) a{4,5}, atab);
        P = cellfun(@(a) a{4,6}, atab);
        F = addnames(array2table(F), vars, vars);
        P = addnames(array2table(P), vars, vars);
        
        EYE = ~eye(size(atab));
        [y,x] = find(EYE);
        ancova_summary = cellfun(@(a) a(2:end, :), atab(EYE), 'uniform', 0);
        ancova_rownames = arrayfun(@(i,j) repmat({[vars{i} '--' vars{j}]}, size(ancova_summary{1}, 1), 1), y, x, 'uniform', 0);
        ancova_summary = vertcat(ancova_summary{:});
        ancova_summary = cell2table(ancova_summary);
        ancova_summary = [vertcat(ancova_rownames{:}) ancova_summary];
        ancova_summary.Properties.VariableNames = {'Combination', 'Summary', 'DF', 'SumSQ', 'MeanSQ', 'F', 'pValue'};
        
        try
            writetable(ancova_summary, [stat_dir 'ancova.xlsx']);
            writetable(F, [stat_dir 'ancova.xlsx'], 'WriteRowNames',true, 'Sheet', 'F');
            writetable(P, [stat_dir 'ancova.xlsx'], 'WriteRowNames',true, 'Sheet', 'P');
%         C = cellfun(@(s) multcompare(s, 'Display', 'off'), stats, 'uniform', 0);
%         p = cellfun(@(c) c(end), C);
        catch ex
            warning('Could not save ancova.xlsx.');
        end
    end
end

function compute_and_save_errorbars(T, output_dir)
    global vars names
    
    T = restrict(T, 3, 3);
    TSun = T(T.sun_shade == 'sonne',:);
    TShade = T(T.sun_shade == 'schatten',:);
    
    [sp, g] = unique(T.sp);
    [~, sorti] = sort(T.HsconRang(g));
    sp = sp(sorti);
    group = T.GroupHscon(g);
    group = group(sorti);
    
    TSunBySp = arrayfun(@(s) TSun(TSun.sp == s,vars), sp, 'uniform', 0);
    TShadeBySp = arrayfun(@(s) TShade(TShade.sp == s,vars), sp, 'uniform', 0);
    
    [~, p] = cellfun(@(sun, shade) ttest2(table2array(sun(:,vars)), ...
                                          table2array(shade(:,vars))), ...
                      TSunBySp, TShadeBySp, 'uniform', 0);
    TMeanSun = cellfun(@(t) nanmean(table2array(t(:,vars)), 1), TSunBySp, 'uniform', 0);
    TMeanShade = cellfun(@(t) nanmean(table2array(t(:,vars)), 1), TShadeBySp, 'uniform', 0);
    TStdSun = cellfun(@(t) nanstd(table2array(t(:,vars)), 1), TSunBySp, 'uniform', 0);
    TStdShade = cellfun(@(t) nanstd(table2array(t(:,vars)), 1), TShadeBySp, 'uniform', 0);
    
    p = cell2mat(p);
    TMeanSun = cell2mat(TMeanSun);
    TMeanShade = cell2mat(TMeanShade);
    TStdSun = cell2mat(TStdSun);
    TStdShade = cell2mat(TStdShade);
                  
    labels = cell(size(p));
    labels(p < 0.05) = {'\textbf{*}'};
    labels(p < 0.01) = {'\textbf{**}'};
    labels(p < 0.001) = {'\textbf{***}'};
    
    [X, Y] = meshgrid(1:numel(vars), 1:numel(sp));
    x_labels = arrayfun(@(x, y) [labels{y,x} char(sp(y))], X, Y, 'uniform', 0);
    y_labels = names;
    
    for i = 1:numel(y_labels)
        errbarplot([TMeanSun(:,i) TMeanShade(:,i)], ...
            [TStdSun(:,i) TStdShade(:,i)], ...
            y_labels{i}, x_labels(:,i), group);
        saveFigure([output_dir vars{i}], [600 400], '-dsvg');
    end
end

function compute_and_save_anovas(T, output_dir)
    global vars
    %%
    T = restrict(T, 3, 3);
    T = T(T.GroupHscon ~= '', :);
    
    models = cellfun(@(trait) [trait ' ~ sp * sun_shade'], vars, 'uniform', 0);
    AV = compute_anova(T, models);
    save_anova(AV, [output_dir 'ANOVA trait - sp x sun_shade.xlsx']);
    
    %%
    groups = categorical(string(T.sp) + string(T.sun_shade));
    meanT = tablemean(T, vars, groups);
    
    models = cellfun(@(trait) [trait ' ~ GroupHscon * sun_shade'], vars, 'uniform', 0);
    AV = compute_anova(meanT, models);
    save_anova(AV, [output_dir 'ANOVA MEAN trait - GroupHscon x sun_shade.xlsx']);
end

function AV = compute_anova(T, model_functions)
    AV = cellfun(@(model) anova(fitlm(T, model)), model_functions, 'uniform', 0);
    
    for i = 1:numel(AV)
        split = strsplit(model_functions{i}, '~');
        trait = strtrim(split{1});
        AV{i}.Properties.RowNames = cellfun(@(name) [trait '-' name], ...
            AV{i}.Properties.RowNames, 'uniform', 0); 
    end

    AV = AV';
    AV = vertcat(AV{:});
end

function save_anova(AV, filename)
    try
        writetable(AV, filename, 'WriteRowNames',true);
    catch
        warning('Could not save anova.xlsx.');
    end
end

function R = restrict(T, minSunSamples, minShadeSamples)
    uniqueSp = unique(T.sp);
    
    U = arrayfun(@(sp) T(T.sp == sp,:), uniqueSp, 'uniform', 0);
    TSun = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'sonne',:), ...
        uniqueSp, 'uniform', 0);
    TShade = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'schatten',:), ...
        uniqueSp, 'uniform', 0);
    
    [countSun,~] = cellfun(@size, TSun);
    [countShade,~] = cellfun(@size, TShade);
    
    sel = countSun >= minSunSamples & countShade >= minShadeSamples;
    R = vertcat(U{sel});
end

function R = tablemean(T, columns, groups)
    [~, g, i] = unique(groups);
    R = T(g, :);
    tmp = arrayfun(@(j) nanmean(table2array(T(j == i, columns)), 1), ...
        (1:numel(g))', 'uniform', 0);
    R(:, columns) = num2cell(vertcat(tmp{:}));
end

function R = tablestd(T, columns, groups)
    [~, g, i] = unique(groups);
    R = T(g, :);
    tmp = arrayfun(@(j) nanstd(table2array(T(j == i, columns)), 1), ...
        (1:numel(g))', 'uniform', 0);
    R(:, columns) = num2cell(vertcat(tmp{:}));
end

function T = addnames(T, names, vars)
    T.Properties.VariableNames = vars;
    T.Properties.RowNames = names;
end

function scatterCompare(X1, Y1, G1, p1, r1, X2, Y2, G2, p2, r2, nameX, nameY, pvalue)
    circle_size = 20;
    font_size = 16;
    ALL = [X1 Y1; X2 Y2];
    MIN = min(ALL);
    MAX = max(ALL);
    border = (MAX - MIN) .* 0.05;
    
    global dark_green light_green

    figure;
    hold on;
    scatter(X1,Y1,circle_size,min(light_green,1),'filled');
    scatter(X2,Y2,circle_size,min(dark_green,1),'filled');
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    h = lsline(ax);
    uistack(h(1), 'up');
    set(h(1),'color',min(dark_green,1));
    if p2 > pvalue
        set(h(1),'linestyle','--');
    end
    set(h(2),'color',min(light_green,1));
    if p1 > pvalue
        set(h(2),'linestyle','--');
    end
    %scatter(MX1,MY1,circle_size,min(color1./255,1),'filled');
    %scatter(MX2,MY2,circle_size,min(color2./255,1),'filled');
    xlim([MIN(1)-border(1) MAX(1)+border(1)]);
    ylim([MIN(2)-border(2) MAX(2)+border(2)]);
    
    xlabel(nameX, 'FontSize', font_size, 'interpreter','latex');
    ylabel(nameY, 'FontSize', font_size, 'interpreter','latex');
    box on;
    
    legend('Sun', 'Shade', ...
           ['r = ' num2str(r1,2) '  p = ' num2str(p1,2)], ...
           ['r = ' num2str(r2,2) '  p = ' num2str(p2,2)]);
    
    hold off;
end

function plotandsavematrix(T, names)
    [~,AX] = plotmatrix(table2array(T));
    vars = T.Properties.VariableNames;
    if nargin < 2
        names = vars;
    end
    [Y,X] = meshgrid(1:numel(vars));
    for xy = [X(:) Y(:)]'
        x = xy(1);
        y = xy(2);
        if x ~= y
            h = lsline(AX(x,y));
            set(h(1),'color','r');
        end
        if y == 1
            xlabel(AX(y,x),names{x}, 'Interpreter', 'latex', 'FontSize', 8);
            set(AX(y,x),'xaxislocation','top');
        end
        if x == numel(vars)
            ylabel(AX(y,x),names{y}, 'Interpreter', 'latex', 'FontSize', 8);
            set(AX(y,x),'yaxislocation','right');
        end
        AX(y,x).TickLabelInterpreter = 'latex';
    end
    
end

function saveFigure(filename, pagesize, format)
    fig = gcf;
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 pagesize];
    print(gcf, filename,format,'-r0');
end

function saveAllSignificantCorr(TSun, TSunSp, TShade, TShadeSp, CorrSun, CorrShade, pvalue, dir)
    [Y,X] = find(triu(ones(size(CorrSun)),1));
    CSun = table2array(CorrSun);
    CShade = table2array(CorrShade);
    for i = [X Y]'
        sun_p = CSun(i(2), i(1));
        sun_r = CSun(i(1), i(2));
        shade_p = CShade(i(2), i(1));
        shade_r = CShade(i(1), i(2));
        if sun_p <= pvalue || shade_p <= pvalue
            nameX = CorrSun.Properties.RowNames{i(1)};
            nameY = CorrSun.Properties.RowNames{i(2)};
            varX = CorrSun.Properties.VariableNames{i(1)};
            varY = CorrSun.Properties.VariableNames{i(2)};
            scatterCompare(...
                TSun(:,i(1)), TSun(:,i(2)), TSunSp, sun_p, sun_r, ...
                TShade(:,i(1)), TShade(:,i(2)), TShadeSp, shade_p, shade_r, ...
                nameX, nameY, pvalue);
            saveFigure([dir varX '-' varY], [400 400], '-dsvg');
        end
    end
end

function errbarplot(mu, sd, figTitle, xLabels, group)
    global light_green dark_green
    
    xticks = 1:numel(xLabels);
    g = unique(group, 'stable');
    groupX = zeros(size(g)); %// central value of each group
    for i = 1:numel(g)
        ticks = xticks(group == g(i)) + (i - 1);
        xticks(group == g(i)) = ticks;
        groupX(i) = mean(ticks);
    end
    
    figure;
    hBar = barwitherr(sd, xticks, mu);
    set(hBar(1), 'FaceColor', light_green);
    set(hBar(2), 'FaceColor', dark_green);
    title(figTitle, 'Interpreter', 'latex', 'FontSize', 14);
    legend('Sun', 'Shade');
    ax = gca;
    ax.XLim = [0.25 max(xticks)+0.75];
    ax.XTick = xticks;
    ax.XTickLabels = xLabels;
    ax.XTickLabelRotation = 45;
    ax.TickLabelInterpreter = 'latex';
    ax.Position = [0.07, 0.16, 0.89, 0.75];
    
    %// Add groups
    yl = ax.YLim;
    groupY = -0.23*(yl(2)-yl(1))+yl(1); %// vertical position of texts. Adjust as needed
    deltaX = -1;
    deltaY = .03; %// controls vertical compression of axis. Adjust as needed
    for i = 1:numel(groupX)
        h = text(groupX(i), groupY, char(g(i)), 'Fontsize', 11, 'Interpreter', 'latex');
        pos = get(h, 'Position');
        ext = get(h, 'Extent');
        pos(1) = pos(1) - ext(3)/2 + deltaX; %// horizontally correct position to make it centered
        set(h, 'Position', pos); %// set corrected position for text
    end
    pos = get(gca, 'position');
    pos(2) = pos(2) + deltaY; %// vertically compress axis to make room for texts
    set(gca, 'Position', pos); %/ set corrected position for axis
    
    mi = min(mu(:));
    ma = max(mu(:));
    d = ma - mi;
    if d / max(abs([mi ma])) < 0.3
        mi = min(mu(:) - sd(:));
        ma = max(mu(:) + sd(:));
        offset = (ma - mi)*0.2;
        ax.YLim = [mi-offset ma+offset];
    end
end
