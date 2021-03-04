function model = convertModelToReversible(model_irrev)
%% model = convertModelToReversible(model_irrev)
% Merges irreversible reaction pairs marked wirh '_r' (created using
% convertModelToIrreversible)
% Input:
%       struct model_irrev:             metabolic model with all
%                                       irreversible reactions
% Output:
%       struct model:                   model with merged irreversible
%                                       reactions


field_names = fieldnames(model_irrev);
rxn_fields = field_names(cellfun(@(x)numel(model_irrev.(x))==numel(model_irrev.rxns),field_names));

remove = [];
% loop over all reactions and find pairs
for i=1:numel(model_irrev.rxns)
    
    % check if the reaction has an associated reverse reaction
    idx = ~cellfun('isempty', regexp(model_irrev.rxns, [model_irrev.rxns{i}, '_r$']));
    if any(idx)
        model_irrev.lb(i) = -1000;
        remove(end+1) = find(idx);
    end
    
end

% remove the reverse reactions and associated field entries
for i=remove'
    for j=1:numel(rxn_fields)
        model_irrev.(rxn_fields{j})(i) = [];
    end
    model_irrev.S(:,i) = [];
end
model = model_irrev;
end