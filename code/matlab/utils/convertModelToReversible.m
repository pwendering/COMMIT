function model = convertModelToReversible(model)
% merges irreversible reaction pairs marked wirh '_r' (created using
% convertModelToIrreversible)

field_names = fieldnames(model);
rxn_fields = field_names(cellfun(@(x)numel(model.(x))==numel(model.rxns),field_names));

remove = [];
% loop over all reactions and find pairs
for i=1:numel(model.rxns)
    
    % check if the reaction has an associated reverse reaction
    idx = ~cellfun('isempty', regexp(model.rxns, [model.rxns{i}, '_r$']));
    if any(idx)
        model.lb(i) = -1000;
        remove(end+1) = find(idx);
    end
    
end

% remove the reverse reactions and associated field entries
for i=remove'
    for j=1:numel(rxn_fields)
        model.(rxn_fields{j})(i) = [];
    end
    model.S(:,i) = [];
end

end