function distMat = svdDistance(models)
% Creates a distance matrix based on the distibution of singular values of
% the stoichiometric matrix scaled to the maximum of each vector of
% singular values. The distance is based on the test statistic of a
% 2-sample Kolmogorov-Smirnov test.
% 
% Input:
%           cell models:            cell array containing stoichiometric
%                                   models with at least field 'S'
% Output:   double distMat:         distance matrix

n = numel(models);

% Check existance of field 'S'
for i=1:n
    if ~isfield(models{i}, 'S')
        error('At least one model does not have a field S containing the stoichiometric matrix')
    end
end

singular_values = cell(n,1);
parfor i=1:n
    singular_values{i} = svd(full(models{i}.S), 'econ');
    %Scaling to the maximum in each set
    singular_values{i} = singular_values{i}/max(singular_values{i});
end

% create the distance matrix
distMat = zeros(n,n);
for i=1:n-1
    for j=i+1:n
        [~,~,stat] = kstest2(singular_values{i}, singular_values{j});
        distMat(i,j) = stat;
        distMat(j,i) = distMat(i,j);
    end
end

end