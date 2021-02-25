function distMat = deadEndsNumberDist(models)
% Creates a distance matrix based on the Jaccard index that compares the
% number of dead-end metabolites that are present in the models. It is important that all
% models have the same namespace to be compared (e.g. KEGG, ModelSEED,
% BiGG,...).
% Input:
%           cell models:            cell array that contains the metabolic
%                                   models, all model structs should at
%                                   least have the 'S', 'mets', 'lb' and 'ub' fields
% Output:   double distMat:         distance matrix

n = numel(models);

for i=1:n
    if ~isfield(models{i}, 'S') || ~isfield(models{i}, 'lb') || ~isfield(models{i}, 'ub') || ~isfield(models{i}, 'mets')
        error('At least one model does not have fields S, mets, lb and/or ub')
    end
end

distMat = zeros(n,n);
for i=1:n-1
    % retrieve the number of dead-ends for model i
    deadEnds_i = numel(detectDeadEnds(models{i}))/numel(models{i}.mets);
    for j=i+1:n
        % retrieve the number of dead-ends for model j
        deadEnds_j = numel(detectDeadEnds(models{j}))/numel(models{j}.mets);
        % Calculation distance
        distMat(i,j) = abs(deadEnds_i-deadEnds_j);
        %distMat(j,i) = distMat(i,j);
    end
end

% Copy the upper right triangle to the lower left triangle to obtain a
% symmetric distance matrix
for i=1:n
    distMat(:,i) = distMat(i,:);
end

end