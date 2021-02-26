function translationDB = loadTranslationDB
options; clearvars -except topDir
% load translation tables
tablesDir = fullfile(topDir, 'data', 'tables');
metTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-met-translation-table.csv');
rxnTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-rxn-translation-table.csv');
metTransTab = readtable(metTransFile, 'ReadVariableNames', true);
rxnTransTab = readtable(rxnTransFile, 'ReadVariableNames', true);
translationDB.metTab = metTransTab;
translationDB.rxnTab = rxnTransTab;
end