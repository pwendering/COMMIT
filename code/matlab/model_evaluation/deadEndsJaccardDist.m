function distMat = deadEndsJaccardDist(models)
% Creates a distance matrix based on the Jaccard index that compares the
% lists of dead-end metabolites that are present in the models. It is important that all
% models have the same namespace to be compared (e.g. KEGG, ModelSEED,
% BiGG,...).
% Input:
%           cell models:            cell array that contains the metabolic
%                                   models, all model structs should at
%                                   least have the 'S', 'lb' and 'ub' fields
% Output:   double distMat:         distance matrix

n = numel(models);

for i=1:n
    if ~isfield(models{i}, 'S') || ~isfield(models{i}, 'lb') || ~isfield(models{i}, 'ub') || ~isfield(models{i}, 'mets')
        error('At least one model does not have fields rxns, mets, lb and/or ub')
    end
end

distMat = zeros(n,n);
for i=1:n-1
    % retrieve the dead-ends for model i
    deadEnds_i = models{i}.mets(detectDeadEnds(models{i}));
    deadEnds_i = regexprep(deadEnds_i, '\[[a-z]\]$', '');
    
    for j=i+1:n
        % retrieve the dead-ends for model j
        deadEnds_j = models{j}.mets(detectDeadEnds(models{j}));
        deadEnds_j = regexprep(deadEnds_j, '\[[a-z]\]$', '');
        % Calculation of the Jaccard distance J(A,B) = intersect(A,B)/union(A,B)
        if ~isempty(deadEnds_i) && ~isempty(deadEnds_j)
            distMat(i,j) = 1-numel(intersect(deadEnds_i, deadEnds_j))/numel(union(deadEnds_i, deadEnds_j));
        else
            distMat(i,j) = 1;
        end

    end
    
end

% Copy the upper right triangle to the lower left triangle to obtain a
% symmetric distance matrix
for i=1:n
    distMat(:,i) = distMat(i,:);
end

end