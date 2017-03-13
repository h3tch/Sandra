function basic

    T = importfile('data.csv');
    T.sp = categorical(T.sp);
    T.sun_shade = categorical(T.sun_shade);
    
    %% SPECIFY OUTPUT FOLDERS
    
    scatter_dir = '../fig/scatter/';
    stat_dir = '../stat/';
    
    %% SELECTE SUB TABLES
    
    uniqueSp = unique(T.sp);
    
    vars = {'Blattdicke_mm', 'Chlorophyll', 'Reissfestigkeit_N', ...
        'DW_FW', 'SLA', 'Stomatadichte', 'd15N14N', 'd13C12C', ...
        'N', 'C', 'PARsat', 'ETR_1500'};
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
    selectedTSun = TSun(sel);
    selectedTShade = TShade(sel);
    
    selTSun = vertcat(selectedTSun{:});
    selTShade = vertcat(selectedTShade{:});
    
    %% SCATTER PLOT COMPARISONS
    
    nameX = 'SLA';
    nameY = 'Blattdicke_mm';
    scatterCompare(...
        table2array(selTSun(:,nameX)), table2array(selTSun(:,nameY)), ...
        table2array(selTShade(:,nameX)), table2array(selTShade(:,nameY)),...
        nameX, nameY);
    saveFigure([scatter_dir nameX '-' nameY], [400 400], '-dsvg');
    
    %% SCATTER PLOT MATRIX
    
    plotandsavematrix(selTSun(:,vars));
    saveFigure([scatter_dir 'scattermatrix_sun'], [1000 1000], '-dpdf');
    plotandsavematrix(selTShade(:,vars));
    saveFigure([scatter_dir 'scattermatrix_shade'], [1000 1000], '-dpdf');
    
    %% PERFORM TTEST
    
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
    
    [h,p] = cellfun(@(u,h) ttest2(table2array(u(:,varSel)),table2array(h(:,varSel))), ...
        selectedTSun, selectedTShade, 'uniform', 0);
    
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
    
    p = cell2mat(p);
    p = array2table(p);
    p.Properties.VariableNames = cellfun(@(i) ['p_',i], vars, 'uniform', 0);
    
    h = cell2mat(h);
    h = array2table(h);
    h.Properties.VariableNames = cellfun(@(i) ['h_',i], vars, 'uniform', 0);
    
    summary = [array2table(selectedSp), mu, sd, num, p, h];
    writetable(summary, [stat_dir 'mean_std_ttest.xlsx']);
    
    %% PERFORM ANOVA per trait
    
    AV = cellfun(@(trait) anova(fitlm(selectedT, [trait ' ~ sp * sun_shade'])), ...
        vars, 'uniform', 0);
    
    for i = 1:numel(AV)
        trait = vars{i};
        AV{i}.Properties.RowNames = cellfun(@(name) [trait '-' name], ...
            AV{i}.Properties.RowNames, 'uniform', 0); 
    end
    
    AV = AV';
    AV = vertcat(AV{:});
    
    AV = [cell2table(AV.Properties.RowNames) AV];
    writetable(AV, [stat_dir 'anova.xlsx']);
    
end

function scatterCompare(X1, Y1, X2, Y2, nameX, nameY)
    circle_size = 13;
    font_size = 10;
    ALL = [X1 Y1; X2 Y2];
    MIN = min(ALL);
    MAX = max(ALL);
    border = (MAX - MIN) .* 0.05;

    figure;
    hold on;
    scatter(X1,Y1,circle_size,[113 183 53]./255,'filled');
    ax = gca;
    lsline(ax);
    scatter(X2,Y2,circle_size,[1 103 52]./255,'filled');
    h = lsline(ax);
    set(h(1),'color','r');
    set(h(2),'color','r');
    xlim([MIN(1)-border(1) MAX(1)+border(1)]);
    ylim([MIN(2)-border(2) MAX(2)+border(2)]);
    
    xlabel(nameX, 'Interpreter', 'none', 'FontSize', font_size);
    ylabel(nameY, 'Interpreter', 'none', 'FontSize', font_size);
    box on;
    
    hold off;
end

function plotandsavematrix(T)
    [~,AX] = plotmatrix(table2array(T));
    vars = T.Properties.VariableNames;
    [Y,X] = meshgrid(1:numel(vars));
    for xy = [X(:) Y(:)]'
        x = xy(1);
        y = xy(2);
        if x ~= y
            h = lsline(AX(x,y));
            set(h(1),'color','r');
        end
        if y == 1
            xlabel(AX(y,x),vars{x}, 'Interpreter', 'none', 'FontSize', 8);
            set(AX(y,x),'xaxislocation','top');
        end
        if x == numel(vars)
            ylabel(AX(y,x),vars{y}, 'Interpreter', 'none', 'FontSize', 8);
            set(AX(y,x),'yaxislocation','right');
        end
    end
    
end

function saveFigure(filename, pagesize, format)
    fig = gcf;
    fig.PaperUnits = 'points';
    fig.PaperPosition = [0 0 pagesize];
    print(gcf, filename,format,'-r0');
end
