% prepare the universal database for conditional gap filling

load('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat')
names = readtable('/stud/wendering/Masterthesis/DATA/tables/MNXref/MET_NAMES_FROM_MNXref.csv',...
    'ReadVariableNames', false);
met_IDs = strtok(dbModel_MNXref_balanced.mets, '[');
idx_mets = cellfun(@(x)find(strcmp(x, names.Var1)), met_IDs,...
    'UniformOutput', false);
names = cellfun(@(x)names.Var2(x), idx_mets, 'UniformOutput', false);
dbModel_MNXref_balanced.metNames = names;

% find membrae-permeable metabolites in the database
dbModel_MNXref_balanced.metPermeability = logical(getPermeabilityWeights(dbModel_MNXref_balanced.mets));
mets_permeable = dbModel_MNXref_balanced.mets(dbModel_MNXref_balanced.metPermeability);
% find reactions in the database that have a transport reaction and their
% metabolites have been predicted to be membrane-permeable
permeablity_rxns = intersect(dbModel_MNXref_balanced.rxns(dbModel_MNXref_balanced.transport),...
    findRxnsFromMets(dbModel_MNXref_balanced, mets_permeable));
% add permeability information to the database model
dbModel_MNXref_balanced.rxnPermeability = zeros(numel(dbModel_MNXref_balanced.rxns),1);
dbModel_MNXref_balanced.rxnPermeability(ismember(dbModel_MNXref_balanced.rxns,permeablity_rxns)) = 1;

save('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat',...
    'dbModel_MNXref_balanced')