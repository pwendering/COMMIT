function distMat = metJaccardDist(models)
% Creates a distance matrix based on the Jaccard index that compares the
% lists of metabolites that are present in the models. It is important that all
% models have the same namespace to be compared (e.g. KEGG, ModelSEED,
% BiGG,...).
% Input:
%           cell models:            cell array that contains the metabolic
%                                   models, all model structs should at least have the 'mets' field
% Output:   double distMat:         distance matrix

n = numel(models);

for i=1:n
    if ~isfield(models{i}, 'mets')
        error('At least one model does not have a field mets containing all metabolites')
    end
end

distMat = zeros(n,n);
for i=1:n-1
    % retrieve the metabolites for model i and remove compartment
    % identifiers
    metabolites_i = models{i}.mets;
    metabolites_i = regexprep(metabolites_i, '\[[a-z]\]$', '');
    
    for j=i+1:n
        % retrieve the metabolites for model j
        metabolites_j = models{j}.mets;
        metabolites_j = regexprep(metabolites_j, '\[[a-z]\]$', '');
        
        % Calculation of the Jaccard distance J(A,B) = intersect(A,B)/union(A,B)
        distMat(i,j) = 1-numel(intersect(metabolites_i, metabolites_j))/numel(union(metabolites_i, metabolites_j));
        distMat(j,i) = distMat(i,j);
    end
end

end