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
    
    startColumn = 17;
    
    nSun = cell2mat(cellfun(@(t) sum(~isnan(table2array(t(:,startColumn:end)))), ...
        selectedTSun, 'uniform', 0));
    nShade = cell2mat(cellfun(@(t) sum(~isnan(table2array(t(:,startColumn:end)))), ...
        selectedTShade, 'uniform', 0));
    
    meanSun = cell2mat(cellfun(@(t) nanmean(table2array(t(:,startColumn:end))), ...
        selectedTSun, 'uniform', 0));
    meanShade = cell2mat(cellfun(@(t) nanmean(table2array(t(:,startColumn:end))), ...
        selectedTShade, 'uniform', 0));
    
    stdSun = cell2mat(cellfun(@(t) nanstd(table2array(t(:,startColumn:end))), ...
        selectedTSun, 'uniform', 0));
    stdShade = cell2mat(cellfun(@(t) nanstd(table2array(t(:,startColumn:end))), ...
        selectedTShade, 'uniform', 0));
    
    httest = cell2mat(cellfun(@(u,h) ttest2(table2array(u(:,startColumn:end)),table2array(h(:,startColumn:end))), ...
        selectedTSun, selectedTShade, 'uniform', 0));
    
    httest = array2table(httest);
    httest.Properties.VariableNames = T.Properties.VariableNames(startColumn:end);
    httest.Properties.RowNames = cellstr(selectedSp);
    
    idx = (1:size(meanSun,2)) * 2 - 1;
    vars = T.Properties.VariableNames(startColumn:end);
    
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
    
    summary = [mu, sd, num];
    
end

