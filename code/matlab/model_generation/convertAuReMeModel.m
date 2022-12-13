function changed_model = convertAuReMeModel(model, translate, translationDB, verbose)
% Convert a metabolic model reconstructed using AuReMe to a form
% that can be used with the FastGapFilling function. For this, field names
% will be changed, added or deleted.
% Input:
%           struct model:               metabolic model that has been
%                                       reconstructed with AuReMe
%           logical translate:          if true, model ids will be translated
%                                       to MNXref ids
%           logical verbose (optional): if true, print warnings and
%                                       progress statements (default: true)
% Output:
%           struct changed_model:       updated AuReMe model

if nargin < 4 || ~islogical(verbose)
    verbose = true;
end

delete = {'metisVersionOfinchiID', 'proteins', 'modelVersion', 'osenseStr',...
    'geneNames'};

% Specify field names that have to be changed
original_names = {'modelID'};
new_names = {'id'};

% Check if all fields that should be changed or deleted are present
check_fields = all(cellfun(@(x)isfield(model, x), [delete, original_names]));
if ~check_fields
    error('The input model does not have all fields that are usually contained in a AuReMe reconstruction')
end

% Change model fields
for i=1:numel(original_names)
   tmp = model.(original_names{i});
   model.(new_names{i}) = tmp;
   model = rmfield(model, original_names{i});
end
clear tmp

% Remove the periplasmal compartment
model.compNames = {'cytosol'; 'extracellular space'};


% Add compartment identifiers to the reactions
to_delete = {};
for i=1:numel(model.rxns)
    if isempty(regexp(model.rxns{i}, '_[ce]$'))
        % get compartment identifiers from all participating metabolites
        idx_mets = logical(model.S(:, i));
        rxn_mets = model.mets(idx_mets);
        comps = extractBetween(rxn_mets, '[', ']');
        
        if numel(unique(comps)) == 1
            % if only one compartment is involved, choose this one
            model.rxns{i} = strcat(model.rxns{i}, '_', comps{1});
        else
            % if it is a transport reaction, chose compartment 'e' for
            % there are should not be  more compartments that 'e' and 'c'
            model.rxns{i} = strcat(model.rxns{i}, '_e');
        end
    end
end

% add 'csense' field to add the type of constraint for each metabolite
model.csense = repmat('E', size(model.S, 1), 1);

% Remove fields from the AuReMe model
model = rmfield(model, delete);

% re-order the field names by length:
field_size = cellfun(@(x)numel(model.(x)), fieldnames(model));
[~,field_order]= sort(field_size);
idx_first = find(ismember(fieldnames(model), {'id', 'description'}));
field_order(ismember(field_order, idx_first)) = [];
field_order = vertcat(idx_first, field_order);
model = orderfields(model, field_order);

if translate
    model = translateModel(model, 'MetaCyc', 'MNXref', translationDB, false, verbose);
end

% unify duplicate metabolite IDs
model = removeDuplicateMets(model);

% unify duplicate reaction IDs
tmp = load(fullfile('data', 'gap-filling', 'database',...
    'Universal-model-MNXref-balanced.mat'));
dbModel = tmp.dbModel_MNXref_balanced;
model = removeDuplicateRxns(model, dbModel);

changed_model = model; 

end
