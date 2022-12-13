function changed_model = convertCarveMeModel(model, translate, translationDB)
% Convert a metabolic model reconstructed using CarveMe to a form
% that can be used with the FastGapFilling function. For this, field names
% will be changed, added or deleted.
% Input:
%           struct model:               metabolic model that has been
%                                       reconstructed with CarveMe
%           logical translate:          if true, model ids will be translated
%                                       to MNXref ids
% Output:
%           struct changed_model:       updated CarveMe model

% Specify fields to delete
delete = {'proteins', 'osenseStr', 'geneNames', 'metFormulas', 'modelVersion',...
    'modelNotes'};

% Specify field names that have to be changed
original_names = {'modelID', 'rxnECNumbers'};
new_names = {'id', 'EC'};

% Check if all fields that should be changed or deleted are present
check_fields = all(cellfun(@(x)isfield(model, x), [delete, original_names]));
if ~check_fields
    error('The input model does not have all fields that are usually contained in a CarveMe reconstruction')
end

% Change model fields
for i=1:numel(original_names)
   tmp = model.(original_names{i});
   model.(new_names{i}) = tmp;
   model = rmfield(model, original_names{i});
end
clear tmp

% Remove 'R_' from the reaction ids
model.rxns = regexprep(model.rxns, '^R_', '');

% change the reaction ids of export reactions by translating the exported
% metabolite ids to the MNXref namespace
ex_rxns = cellfun(@(x)contains(x, 'EX_'), model.rxns);
tmp_suffix = cellfun(@(x)regexp(x, '_.$', 'match'), model.rxns(ex_rxns));
tmp_mets = extractBetween(model.rxns(ex_rxns), 'EX_', '_e');
tmp_mets = translateIDs(tmp_mets, 'met', translationDB.metTab, 'BiGG', 'MNXref');
model.rxns(ex_rxns) = strcat('EX_', tmp_mets, tmp_suffix);
clear ex_rxns tmp_suffix tmp_mets

% change the reaction ids of sink reactions by translating the exported
% metabolite ids to the MNXref namespace
sink_rxns = cellfun(@(x)contains(x, 'sink_'), model.rxns);
tmp_suffix = cellfun(@(x)regexp(x, '_.$', 'match'), model.rxns(sink_rxns));
tmp_mets = extractBetween(model.rxns(sink_rxns), 'sink_', '_c');
tmp_mets = translateIDs(tmp_mets, 'met',  translationDB.metTab, 'BiGG', 'MNXref');
model.rxns(sink_rxns) = strcat('EX_', tmp_mets, tmp_suffix);
clear sink_rxns tmp_suffix tmp_mets

% change the compartment id
model.mets = regexprep(model.mets, '\[C_', '\[');

% Remove the periplasmal compartment
model.comps = {'c'; 'e'};
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
            if any(contains(comps, 'p'))
                comps = {'e'};
            end
            % if only one compartment is involved, choose this one
            model.rxns{i} = strcat(model.rxns{i}, '_', comps{1});
        else
            % if it is a transport reaction, chose compartment 'e' for
            % there are should not be  more compartments that 'e' and 'c'
            if all(contains({'p', 'e'}, comps))
                % if it is a transport reaction between 'e' and 'p', remove
                % the reaction
                to_delete = vertcat(to_delete, model.rxns{i});
            else
                model.rxns{i} = strcat(model.rxns{i}, '_e');
            end
        end
    end
end

% move all metabolites remaining in the periplasm to the extra-cellular
% compartment (change in every reaction and delete the old ids)
old_met_idx = find(contains(model.mets, '[p]'));
old_mets = model.mets(contains(model.mets, '[p]'));
new_mets = regexprep(old_mets, '\[p\]$', '\[e\]');
for i=1:numel(old_mets)
    rxns_met = findRxnsFromMets(model, old_mets{i});
    for j=1:numel(rxns_met)
        if ~ismember(new_mets{i}, model.mets)
            name = model.metNames{old_met_idx(i)};
            formula = model.metFormulas(old_met_idx(i));
            note = model.metNotes(old_met_idx(i));
            model = addMetabolite(model, new_mets{i}, 'metName', name);
            model.metFormulas(end) = formula;
            model.metNotes(end) = note;
        end
        model = changeRxnMets(model, old_mets{i}, new_mets{i}, rxns_met{j});
    end
end
model = removeMetabolites(model, old_mets);

% remove transport reactions from the extra-cellular space to the cytosol
model = removeRxns(model, to_delete);

% add 'csense' field to add the type of constraint for each metabolite
model.csense = repmat('E', size(model.S, 1), 1);

% Change the id of the biomass reaction
idx = contains(model.rxns, 'Growth');
model.rxns(idx) = {'BIOMASS_Reaction'};

% Set the objective to Biomass
model.c(idx) = 1;

% change '_' to '.' in the gene names
model.genes = cellfun(@(x)regexprep(x, '_', '\.'), model.genes,...
    'UniformOutput', false);
model.genes = cellfun(@(x)regexprep(x, '\.', '\|', 'once'), model.genes,...
            'UniformOutput', false);
% Remove fields from the CarveMe model
model = rmfield(model, delete);

% re-order the field names by size:
field_size = cellfun(@(x)numel(model.(x)), fieldnames(model));
[~,field_order]= sort(field_size);
idx_first = find(ismember(fieldnames(model), {'id', 'description'}));
field_order(ismember(field_order, idx_first)) = [];
field_order = vertcat(idx_first, field_order);
model = orderfields(model, field_order);

if translate
    model = translateModel(model, 'BiGG', 'MNXref', translationDB);
    clear translationDB
    % reactions from rxnNotes
    for i=1:numel(model.rxns)
        if ~contains(model.rxns{i}, 'MNXR')
            new_rxn = regexp(model.rxnNotes{i}, 'MNXR[0-9]*', 'match');
            if ~isempty(new_rxn)
                tmp_suffix = regexp(model.rxns{i}, '_.$', 'match');
                model.rxns(i) = strcat(new_rxn, tmp_suffix);
            end
        end
    end
    fprintf('%3.2f%% of reactions are now in MNXref namespace',...
        100*sum(contains(model.rxns, 'MNXR'))/numel(model.rxns));
    clear new_rxn tmp_suffix
    
    % Metabolites from metNotes
    for i=1:numel(model.mets)
        if ~contains(model.mets{i}, 'MNXM')
            new_met = regexp(model.metNotes{i}, 'MNXM[0-9]*', 'match');
            if ~isempty(new_met)
                comp =  regexp(model.mets{i}, '\[.\]', 'match');
                model.rxns(i) = strcat(new_met, comp);
            end
        end
    end
    clear new_met comp
    fprintf('\n%3.2f%% of metabolites are now in MNXref namespace\n',...
        100*sum(contains(model.mets, 'MNXM'))/numel(model.mets));
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