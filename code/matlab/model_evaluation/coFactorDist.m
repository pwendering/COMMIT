function distMat = coFactorDist(models, cofactorList)
% Creates a distance matrix based on the usage of cofactors provided as a
% list in the same namespace as the model is in. 
% The abundances are scaled to the number of reactions in the respective model.
% The distance is calculated with spearman correlation.
% Input:
%               cell models:            array that contains the metabolic
%                                       models that a least contain the field 'rxns'
%               cell cofactorList:      array that contains the cofactors
%                                       the are in the same namespace as 
%                                       the metabolites in the model
% Output:
%               double distMar:     	distance matrix
%    
n = numel(models);
m = numel(cofactorList);

for i=1:n
    if ~isfield(models{i}, 'rxns')
        error('At least one model does not have a field rxns containing Enzyme Commission numbers')
    end
end

% Get a matrix for the usage of cofactors (abundance)
cfUsage = cofactorUsage(models, cofactorList);

% scale each row to the number of reactions in the model
for i=1:n
    cfUsage(i,:) = cfUsage(i,:) / numel(models{i}.rxns);
end

distMat = zeros(n,n);

for i=1:n-1
    % Get a vector of abundances for cofactors
    v_i = cfUsage(i,:)';
    for j=i+1:n
        % Get a vector of abundances for cofactors
        v_j = cfUsage(j,:)';
        try
            distMat(i,j) = 1 - corr(v_i,v_j,'type','Spearman');
            distMat(j,i) = distMat(i,j);
        catch
                disp(v_i)
                disp(v_j)
                error('corr not possible')
        end
    end
end

end