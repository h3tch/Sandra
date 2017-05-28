function basic

    T = importfile('data.csv');
    T.sp = categorical(T.sp);
    T.sun_shade = categorical(T.sun_shade);
    
    %% SPECIFY OUTPUT FOLDERS
    
    barplot_dir = 'fig/bar/';
    scatter_dir = 'fig/scatter/';
    stat_dir = 'stat/';
    mkdir(stat_dir);
    mkdir(barplot_dir);
    mkdir(scatter_dir);
    global dark_green light_green
    dark_green = [113 183 53]./255;
    light_green = [1 103 52]./255;
    
    %% SELECTE SUB TABLES
    
    [uniqueSp, uniqueIdx] = unique(T.sp);
    
    vars = {'Blattdicke_mm', 'Chlorophyll', 'Reissfestigkeit_N', ...
        'DW_FW', 'SLA', 'Stomatadichte', 'd15N14N', 'd13C12C', ...
        'N', 'C', 'PARsat', 'ETR_1500', 'Hzuwachs'};
    names = {'Blattdicke ($mm$)', 'Chlorophyll (SPAD)', 'Reissfestigkeit ($N$)', ...
        'LDMC ($mg/g$)', 'SLA  ($mm^2/mg$)', 'Stomatadichte ($n/mm^2$)', ...
        '$d~15N/14N$', '$d~13C/12C$', 'N ($\%$)', 'C ($\%$)', 'PARsat', ...
        'ETR1500 ($\mu$mol $m^{-2}~s^{-1}$)', 'H\"ohenzuwachs ($cm/year$)'};
    varSel = ismember(T.Properties.VariableNames, vars);
    
    U = arrayfun(@(sp) T(T.sp == sp,:), uniqueSp, 'uniform', 0);
    TSun = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'sonne',:), ...
        uniqueSp, 'uniform', 0);
    TShade = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'schatten',:), ...
        uniqueSp, 'uniform', 0);
    
    [countSun,~] = cellfun(@size, TSun);
    [countShade,~] = cellfun(@size, TShade);
    
    sel = countSun > 3 & countShade > 3;
    selectedT = vertcat(U{sel});
    selectedSp = uniqueSp(sel);
    selectedSpGroup = T.Group(uniqueIdx);
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
    
    if 0
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
    
    if 0
%         TSunSp = selTSun(:,'sp');
%         TShadeSp = selTShade(:,'sp');
        
        pvalue = 0.05;
        saveAllSignificantCorr(TSun, [], TShade, [], CorrSun, CorrShade, pvalue, scatter_dir);
        close all;
    end
    
    %% SCATTER PLOT MATRIX
    
    if 0
        plotandsavematrix(selTSun(:,vars), names);
        saveFigure([scatter_dir 'scattermatrix_sun'], [1000 1000], '-dpdf');
        plotandsavematrix(selTShade(:,vars), names);
        saveFigure([scatter_dir 'scattermatrix_shade'], [1000 1000], '-dpdf');
        
        close all
    end
    
    %% PERFORM TTEST
    
    if 0
        nSun = cell2mat(cellfun(@(t) sum(~isnan(table2array(t(:,varSel)))), ...
            selectedTSun, 'uniform', 0));
        nShade = cell2mat(cellfun(@(t) sum(~isnan(table2array(t(:,varSel)))), ...
            selectedTShade, 'uniform', 0));

        meanSun = cell2mat(cellfun(@(t) nanmean(table2array(t(:,varSel))), ...
            selectedTSun, 'uniform', 0));
        meanShade = cell2mat(cellfun(@(t) nanmean(table2array(t(:,varSel))), ...
            selectedTShade, 'uniform', 0));

        stdSun = cell2mat(cellfun(@(t) nanstd(table2array(t(:,varSel))), ...
            selectedTSun, 'uniform', 0));
        stdShade = cell2mat(cellfun(@(t) nanstd(table2array(t(:,varSel))), ...
            selectedTShade, 'uniform', 0));
        
        % sort by group
        [~,sortIdx] = sort(selectedSpGroup);
        selectedSp = selectedSp(sortIdx);
        meanSun = meanSun(sortIdx,:);
        meanShade = meanShade(sortIdx,:);
        stdSun = stdSun(sortIdx,:);
        stdShade = stdShade(sortIdx,:);
        nSun = nSun(sortIdx,:);
        nShade = nShade(sortIdx,:);

        [h,p] = cellfun(@(u,h) ttest2(table2array(u(:,varSel)),table2array(h(:,varSel))), ...
            selectedTSun(sortIdx), selectedTShade(sortIdx), 'uniform', 0);
        p = cell2mat(p);
        labels = cell(size(p));
        labels(p < 0.05) = {'\textbf{*}'};
        labels(p < 0.01) = {'\textbf{**}'};
        labels(p < 0.001) = {'\textbf{***}'};
        
        for i = 1:size(meanSun,2)
            errbarplot([meanSun(:,i) meanShade(:,i)], ...
                [stdSun(:,i) stdShade(:,i)], ...
                names{i}, arrayfun(@(sp,l) [l{1} '\textsf{' char(sp) '}'], selectedSp, labels(:,i), 'uniform', 0));
            saveFigure([barplot_dir vars{i}], [600 400], '-dsvg');
        end

        idx = (1:size(meanSun,2)) * 2 - 1;

        mu = zeros(size(meanSun) .* [1 2]);
        mu(:,idx) = meanSun;
        mu(:,idx + 1) = meanShade;
        mu = array2table(mu);
        mu.Properties.VariableNames(idx) = cellfun(@(i) ['mu_sun_',i], vars, 'uniform', 0);
        mu.Properties.VariableNames(idx + 1) = cellfun(@(i) ['mu_shade_',i], vars, 'uniform', 0);

        sd = zeros(size(stdSun) .* [1 2]);
        sd(:,idx) = stdSun;
        sd(:,idx + 1) = stdShade;
        sd = array2table(sd);
        sd.Properties.VariableNames(idx) = cellfun(@(i) ['sd_sun_',i], vars, 'uniform', 0);
        sd.Properties.VariableNames(idx + 1) = cellfun(@(i) ['sd_shade_',i], vars, 'uniform', 0);

        num = zeros(size(stdSun) .* [1 2]);
        num(:,idx) = nSun;
        num(:,idx + 1) = nShade;
        num = array2table(num);
        num.Properties.VariableNames(idx) = cellfun(@(i) ['num_sun_',i], vars, 'uniform', 0);
        num.Properties.VariableNames(idx + 1) = cellfun(@(i) ['num_shade_',i], vars, 'uniform', 0);

        p = array2table(p);
        p.Properties.VariableNames = cellfun(@(i) ['p_',i], vars, 'uniform', 0);

        h = cell2mat(h);
        h = array2table(h);
        h.Properties.VariableNames = cellfun(@(i) ['h_',i], vars, 'uniform', 0);

        summary = [array2table(selectedSp), mu, sd, num, p, h];
        writetable(summary, [stat_dir 'mean_std_ttest.xlsx'], 'WriteRowNames',true);
        
        close all
    end
    
    %% PERFORM ANOVA per trait
    
    if 1
        AV = cellfun(@(trait) anova(fitlm(selectedT, [trait ' ~ sp * sun_shade'])), ...
            vars, 'uniform', 0);

        for i = 1:numel(AV)
            trait = vars{i};
            AV{i}.Properties.RowNames = cellfun(@(name) [trait '-' name], ...
                AV{i}.Properties.RowNames, 'uniform', 0); 
        end

        AV = AV';
        AV = vertcat(AV{:});

        %AV = [cell2table(AV.Properties.RowNames) AV];
        try
            writetable(AV, [stat_dir 'anova.xlsx'], 'WriteRowNames',true);
        catch
            warning('Could not save anova.xlsx.');
        end
    end
    
    %% PERFORM ANCOVA
    
    if 1
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

function B = sortmat(A, sortIdx)
    B = A(:,sortIdx);
    B = B(sortIdx,:);
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
    scatter(X1,Y1,circle_size,min(dark_green,1),'filled');
    scatter(X2,Y2,circle_size,min(light_green,1),'filled');
    ax = gca;
    ax.TickLabelInterpreter = 'latex';
    h = lsline(ax);
    set(h(1),'color',min(light_green,1));
    if p1 > pvalue
        set(h(1),'linestyle','--');
    end
    set(h(2),'color',min(dark_green,1));
    if p2 > pvalue
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
        nameX = CorrSun.Properties.RowNames{i(1)};
        nameY = CorrSun.Properties.RowNames{i(2)};
        varX = CorrSun.Properties.VariableNames{i(1)};
        varY = CorrSun.Properties.VariableNames{i(2)};
        scatterCompare(...
            TSun(:,i(1)), TSun(:,i(2)), TSunSp, CSun(i(2), i(1)), CSun(i(1),i(2)), ...
            TShade(:,i(1)), TShade(:,i(2)), TShadeSp, CShade(i(2), i(1)), CShade(i(1),i(2)), ...
            nameX, nameY, pvalue);
        saveFigure([dir varX '-' varY], [400 400], '-dsvg');
    end
end

function errbarplot(mu, sd, figTitle, xLabels, groupLabels)
    global light_green dark_green
    figure;
    hBar = barwitherr(sd, mu);
    set(hBar(2), 'FaceColor', light_green);
    set(hBar(1), 'FaceColor', dark_green);
    title(figTitle, 'Interpreter', 'latex', 'FontSize', 14);
    legend('Sun', 'Shade');
    ax = gca;
    ax.XLim = [0.25 numel(xLabels)+0.75];
    ax.XTick = 1:numel(xLabels);
    ax.XTickLabels = xLabels;
    ax.XTickLabelRotation = 45;
    ax.TickLabelInterpreter = 'latex';
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
