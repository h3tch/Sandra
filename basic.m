function basic

    T = importfile('data_2017.03.01.csv');
    T.sp = categorical(T.sp);
    T.sun_shade = categorical(T.sun_shade);
    
    uniqueSp = unique(T.sp);
    
    TSun = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'sonne',:), ...
        uniqueSp, 'uniform', 0);
    TShade = arrayfun(@(sp) T(T.sp == sp & T.sun_shade == 'schatten',:), ...
        uniqueSp, 'uniform', 0);
    
    [countSun,~] = cellfun(@size, TSun);
    [countShade,~] = cellfun(@size, TShade);
    
    sel = countSun > 3 & countShade > 3;
    selectedSp = uniqueSp(sel);
    selectedTSun = TSun(sel);
    selectedTShade = TShade(sel);
    
    meanSun = cell2mat(cellfun(@(t) nanmean(table2array(t(:,6:end))), ...
        selectedTSun, 'uniform', 0));
    meanShade = cell2mat(cellfun(@(t) nanmean(table2array(t(:,6:end))), ...
        selectedTShade, 'uniform', 0));
    
    stdSun = cell2mat(cellfun(@(t) nanstd(table2array(t(:,6:end))), ...
        selectedTSun, 'uniform', 0));
    stdShade = cell2mat(cellfun(@(t) nanstd(table2array(t(:,6:end))), ...
        selectedTShade, 'uniform', 0));
    
    httest = cell2mat(cellfun(@(u,h) ttest2(table2array(u(:,6:end)),table2array(h(:,6:end))), ...
        selectedTSun, selectedTShade, 'uniform', 0));
    
    httest = array2table(httest);
    httest.Properties.VariableNames = T.Properties.VariableNames(6:end);
    httest.Properties.RowNames = cellstr(selectedSp);
    
end

