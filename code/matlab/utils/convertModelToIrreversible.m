function model_irrev = convertModelToIrreversible(model)
%% model_irrev = convertModelToIrreversible(model)
% Splits all reversible reactions in a metabolic model into two
% irreversible reaction and updates all associated fields
% Input:
%       struct model:           metabolic model
% Output:
%       struct model_irrev:     all-irreversible model


% get indices of all reversible reactions (both lb and ub != 0)
rev = find(model.lb<0&model.ub>0);

% generate a new matrix containing the new unidirectional reactions
irrev_matrix = -model.S(:, rev);

% add the new matrix to the original model
model.S = [model.S, irrev_matrix];

% change the upper and lower boundaries for the now split reversible reactions
model.lb(rev) = 0;
model.ub(rev) = 1000;

% add boundaries for the added reactions
model.lb = vertcat(model.lb, zeros(numel(rev), 1));
model.ub = vertcat(model.ub, repmat(1000, numel(rev), 1));

% construct new ids by adding '_r' to the new irreversible reactions
irrev_ids = strcat(model.rxns(rev), '_r');

% add the new ids to the 'rxns' field
model.rxns = vertcat(model.rxns, irrev_ids);

% update reaction-related fields
if isfield(model, 'transport')
    ids_transport = model.rxns(logical(model.transport));
    model.transport = contains(model.rxns, ids_transport);
end

if isfield(model, 'c')
    model.c = vertcat(model.c, model.c(rev));
end

try
    if isfield(model, 'EC')
        model.EC = vertcat(model.EC, model.EC(rev));
    end
catch
    model = rmfield(model, 'EC');
end

if isfield(model, 'rxnNames')
    model.rxnNames = vertcat(model.rxnNames, model.rxnNames(rev));
end

if isfield(model, 'rxnNotes')
    model.rxnNotes = vertcat(model.rxnNotes, model.rxnNotes(rev));
end

if isfield(model, 'rules')
    model.rules = vertcat(model.rules, model.rules(rev));
end

if isfield(model, 'subSystems')
    model.subSystems = vertcat(model.subSystems, model.subSystems(rev));
end

if isfield(model, 'scores')
    model.scores = vertcat(model.scores, model.scores(rev));
end

if isfield(model, 'rxnPermeability')
    ids_permeability = model.rxns(logical(model.rxnPermeability));
    model.rxnPermeability = contains(model.rxns, ids_permeability);
end
model_irrev = model;

end
