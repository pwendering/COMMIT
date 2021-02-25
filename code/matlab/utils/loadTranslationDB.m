function translationDB = loadTranslationDB
% load translation tables
tablesDir = '/stud/wendering/Masterthesis/DATA/tables';
metTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-met-translation-table.csv');
rxnTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-rxn-translation-table.csv');
metTransTab = readtable(metTransFile, 'ReadVariableNames', true);
rxnTransTab = readtable(rxnTransFile, 'ReadVariableNames', true);
translationDB.metTab = metTransTab;
translationDB.rxnTab = rxnTransTab;
clear metTransTab rxnTransTab tablesDir metTransFile rxnTransFile
end