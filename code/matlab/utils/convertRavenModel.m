function changed_model = convertRavenModel(model, translate, translationDB)
% Convert a metabolic model reconstructed using the RAVEN toolbox to a form
% that can be used with the FastGapFilling function. For this, field names
% will be changed, added or deleted.
% Input:
%           struct model:           metabolic model that has been
%                                   reconstructed with RAVEN
%           logical translate:      if true, model ids will be translated
%                                   to MNXref ids
% Output:
%           struct changed_model:   updated RAVEN model

% Specify fields to delete
delete = {'rxnMiriams', 'rxnReferences', 'rev', 'equations',...
    'inchis', 'metMiriams', 'metComps', 'rxnFrom', 'geneFrom', 'grRulesKEGG',...
    'metFrom', 'metFormulas', 'rxnConfidenceScores', 'rxnGeneMat'};

% Specify field names that have to be changed
original_names = {'eccodes', 'grRules' };
new_names = {'EC', 'rules'};

% Check if all fields that should be changed or deleted are present
check_fields = all(cellfun(@(x)isfield(model, x), [delete, original_names]));
if ~check_fields
    error('The input model does not have all fields that are usually contained in a RAVEN reconstruction')
end

% Change model fields
for i=1:numel(original_names)
   tmp = model.(original_names{i});
   model.(new_names{i}) = tmp;
   model = rmfield(model, original_names{i});
end
clear tmp

% Combine fields
model.rxnNotes = strcat('Reference:', {' '}, model.rxnReferences,...
    '; From:', {' '}, model.rxnFrom);

model.metNotes = strcat('From:', {' '}, model.metFrom,...
    '; InChI:', {' '}, model.inchis);

% Change model Compartment from 's'/ 'System' to 'c0' 'Cytosol'
if numel(model.comps) < 3
    % add the extracellular compartment in any case
    model.comps = {'c', 'e'};
    model.compNames = {'Cytosol', 'Extracellular'};
    for i=1:numel(model.mets)
        if model.metComps(i) == 1
            model.mets{i} = strcat(model.mets{i}, '[', model.comps{1}, ']');
        else
            model.mets{i} = strcat(model.mets{i}, '[', model.comps{2}, ']');
        end
    end
elseif numel(model.comps) == 3
    model.comps = {'c', 'e', 'p'};
    model.compNames = {'Cytosol', 'Extracellular', 'Periplasm'};
    for i=1:numel(model.mets)
        if model.metComps(i) == 1
            model.mets{i} = strcat(model.mets{i}, '[', model.comps{1}, ']');
        elseif model.metComps(i) == 2
            model.mets{i} = strcat(model.mets{i}, '[', model.comps{2}, ']');
        else
            model.mets{i} = strcat(model.mets{i}, '[', model.comps{3}, ']');
        end
    end
else 
    error('Too many compartment identifiers for a bacterial model')
end
% There are no specified compartments reactions so all of them are assigned
% to the cytosol
model.rxns = strcat(model.rxns, '_c');

% Change the format of E.C. numbers
model.EC = regexprep(model.EC, 'ec-code/', '');
model.EC = regexprep(model.EC, ';', '|');

% add 'csense' field to add the type of constraint for each metabolite
model.csense = repmat('E', size(model.S, 1), 1);

% Remove fields from the RAVEN model
model = rmfield(model, delete);

% re-order the field names by length:
field_size = cellfun(@(x)numel(model.(x)), fieldnames(model));
[~,field_order]= sort(field_size);
idx_first = find(ismember(fieldnames(model), {'id', 'description'}));
field_order(ismember(field_order, idx_first)) = [];
field_order = vertcat(idx_first, field_order);
model = orderfields(model, field_order);

if translate
    % first translate KEGG IDs, then MetaCyc IDs
    model = translateModel(model, 'MetaCyc', 'MNXref', translationDB);
    model = translateModel(model, 'KEGG', 'MNXref', translationDB);
end

% unify duplicate metabolite IDs
model = removeDuplicateMets(model);

% unify duplicate reaction IDs
model = removeDuplicateRxns(model);

changed_model = model;
end