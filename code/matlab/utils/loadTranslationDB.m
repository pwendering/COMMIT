function translationDB = loadTranslationDB
%% translationDB = loadTranslationDB
% load translation tables for metabolites and reactions
% Output:
%   struct translationDB:           structure with fields 'metTab' and
%                                   'rxnTab'

tablesDir = fullfile('data', 'tables');
metTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-met-translation-table.csv');
rxnTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-rxn-translation-table.csv');
metTransTab = readtable(metTransFile, 'ReadVariableNames', true,...
    'Format', '%s\t%s\t%s\t%s\t%s\t%s\t%s');
rxnTransTab = readtable(rxnTransFile, 'ReadVariableNames', true,...
    'Format', '%s\t%s\t%s\t%s\t%s\t%s\t%s');
translationDB.metTab = metTransTab;
translationDB.rxnTab = rxnTransTab;
end