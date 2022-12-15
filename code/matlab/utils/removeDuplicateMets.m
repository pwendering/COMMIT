function changed_model = removeDuplicateMets(model)
%% changed_model = removeDuplicateMets(model)
% Combines rows of metabolites in the stoichiometric matrix, which share 
% the same identifier
% Input:
%           struct model:           metabolic model
% Output:   
%           struct changed_model:   updated model

% first occurrences of metabolites
[~, ia, ~] = unique(model.mets, 'stable');

% indices of duplicate metabolites
to_delete = model.mets(setdiff(1:numel(model.mets), ia));
remove = {''};

for i=1:numel(to_delete)
    % indices of the repeated occurrences of the same metabolite
    idx = find(contains(model.mets, to_delete(i)));
    % add the coefficients of the other metabolites (in other reactions)
    % to this id
    model.S(idx(1), :) = sum(model.S(idx, :), 1);
    for j=2:numel(idx)
        % rename the id to be able to remove it savely
        model.mets(idx(j)) = strcat(model.mets(idx(j)), '_', num2str(j));
        % remove the metabolite from the model
        remove = vertcat(remove, model.mets(idx(j)));
    end
end

changed_model = removeMetabolites(model, remove, 0);

end