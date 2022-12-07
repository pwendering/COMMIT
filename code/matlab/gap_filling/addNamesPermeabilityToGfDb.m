function db_out = addNamesPermeabilityToGfDb(db_in, mnxref_dir)
%% db_out = addNamesPermeabilityToGfDb(db_in, mnxref_dir)
% Add metabolite names and permeability information for transport
% reactions to the gap-filling database.
% Input:
%   struct db_in:           input gap-filling database
%   char mnxref_dir:        path directory where MetaNetX files are located
% Output:
%   struct db_out:          gap-filling database with metabolite names and
%                           permeability

if ~isfile(fullfile(mnxref_dir, 'MET_NAMES_FROM_MNXref.csv'))
    error('Could not read file: %s',...
        fullfile(mnxref_dir, 'MET_NAMES_FROM_MNXref.csv'))
end

db_out = db_in;

% read metabolite names
met_names = readtable(...
    fullfile(mnxref_dir, 'MET_NAMES_FROM_MNXref.csv'), ...
    'ReadVariableNames', false);

% remove compartment information from metabolite IDs
met_ids = strtok(db_in.mets, '[');

% find matching metabolite names
idx_mets = cellfun(@(x)find(strcmp(x, met_names.Var1)), met_ids,...
    'UniformOutput', false);
met_names = cellfun(@(x)met_names.Var2(x), idx_mets, 'UniformOutput', false);
db_out.metNames = met_names;

% find membrane-permeable metabolites in the database
db_out.metPermeability = logical(getPermeabilityWeights(db_in.mets));
mets_permeable = db_out.mets(db_out.metPermeability);

% find reactions in the database that have a transport reaction and their
% metabolites have been predicted to be membrane-permeable
permeablity_rxns = intersect(db_in.rxns(db_in.transport),...
    findRxnsFromMets(db_in, mets_permeable));

% add permeability information to the database model
db_out.rxnPermeability = zeros(numel(db_in.rxns),1);
db_out.rxnPermeability(ismember(db_in.rxns,permeablity_rxns)) = 1;

end