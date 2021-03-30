function changed_model = convertKBaseModel(model, translate, translationDB)
% Convert a metabolic model reconstructed using the KBase web service to a form
% that can be used with the FastGapFilling function. For this, field names
% will be changed, added or deleted.
% Input:
%           struct model:               metabolic model that has been
%                                       reconstructed with KBase
%           logical translate:          if true, model ids will be translated
%                                       to MNXref ids
% Output:
%           struct changed_model:       updated KBase model

% Specify fields to delete
delete = {'proteins', 'osenseStr'};

% Specify field names that have to be changed
original_names = {'modelID', 'modelName' };
new_names = {'id', 'description'};

% Check if all fields that should be changed or deleted are present
check_fields = all(cellfun(@(x)isfield(model, x), [delete, original_names]));
if ~check_fields
    error('The input model does not have all fields that are usually contained in a KBase reconstruction')
end

% Change model fields
for i=1:numel(original_names)
   tmp = model.(original_names{i});
   model.(new_names{i}) = tmp;
   model = rmfield(model, original_names{i});
end
clear tmp

% Remove '0' from the model compartments
model.comps = regexprep(model.comps, '0', '');
model.compNames = regexprep(model.compNames, '_0', '');

% Remove '0' from the reaction ids, reaction names and metabolite ids
model.rxns = regexprep(model.rxns, '0$', '');
model.rxns = regexprep(model.rxns, '^R_', '');

% change the reaction ids of export reactions by translating the exported
% metabolite ids to the MNXref namespace
ex_rxns = cellfun(@(x)contains(x, 'EX_'), model.rxns);
tmp_suffix = cellfun(@(x)regexp(x, '_.$', 'match'), model.rxns(ex_rxns));
tmp_mets = regexp(model.rxns(ex_rxns), 'cpd[0-9]*', 'match');
tmp_mets = translateIDs(tmp_mets, 'met', translationDB.metTab, 'ModelSEED', 'MNXref');
model.rxns(ex_rxns) = strcat('EX_', tmp_mets, tmp_suffix);

model.rxnNames = regexprep(model.rxnNames, '0$', '');
model.mets = regexprep(model.mets, '0\]', '\]');
clear ex_rxns tmp_suffix tmp_mets

% Change the id of the biomass reaction
model.rxns{logical(model.c)} = 'BIOMASS_Reaction';

% Remove fields from the KBase model
changed_model = rmfield(model, delete);

% re-order the field names by length:
field_size = cellfun(@(x)numel(changed_model.(x)), fieldnames(changed_model));
[~,field_order]= sort(field_size);
idx_first = find(ismember(fieldnames(changed_model), {'id', 'description'}));
field_order(ismember(field_order, idx_first)) = [];
field_order = vertcat(idx_first, field_order);
changed_model = orderfields(changed_model, field_order);

if translate
    changed_model = translateModel(changed_model, 'ModelSEED', 'MNXref', translationDB);
end


end