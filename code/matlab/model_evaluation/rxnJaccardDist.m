function distMat = rxnJaccardDist(models)
% Creates a distance matrix based on the Jaccard index that compares the
% lists of reactions that are present in the models. It is important that all
% models have the same namespace to be compared (e.g. KEGG, ModelSEED,
% BiGG, MNXref,...).
% Input:
%           cell models:            cell array that contains the metabolic
%                                   models, all model structs should at least have the 'rxns' field
% Output:   double distMat:         distance matrix

n = numel(models);

for i=1:n
    if ~isfield(models{i}, 'rxns')
        error('At least one model does not have a field rxns containing all reactions')
    end
end

distMat = zeros(n);
for i=1:n-1
    % retrieve the reactions for model i
    reactions_i = models{i}.rxns;
    
    % Exclude exchange reactions and remove compartment identifier
    exclude = contains(reactions_i, 'EX_');
    reactions_i = reactions_i(~exclude);
    reactions_i = regexprep(reactions_i, '_[a-z]$', '');

    for j=i+1:n
        % retrieve the reactions for model j
        reactions_j = models{j}.rxns;
        exclude = contains(reactions_j, 'EX_');
        reactions_j = reactions_j(~exclude);
        reactions_j = regexprep(reactions_j, '_[a-z]$', '');
        
        % Calculation of the Jaccard distance J(A,B) = |intersect(A,B)| / |union(A,B)|
        distMat(i,j) = 1-numel(intersect(reactions_i, reactions_j))/numel(union(reactions_i, reactions_j));
        distMat(j,i) = distMat(i,j);
    end
end

end