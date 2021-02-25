function distMat = N1RvDist(models, blackList, biomass, type)
% Create a distance matrix that compares the connectivity of reactions or
% Enzyme Commission (EC) numbers based by calculating the Rv-coefficient.
% Input:
%           cell models:            array containing the metabolic models
%                                   that at least contain the field 'rxns', 
%                                   'mets','EC' and 'S'
%           cell blacklist:         metabolite ids that should be excluded
%           char type:              type of graph that should be
%                                   constructed ('ec' or 'rxns'(default))
%           char biomass:           reaction id of the biomass reaction
% Output:   double distMat:         distance matrix

n = numel(models);

if ~exist('type', 'var') || ~any(ismember(type, {'EC' 'rxns'}))
    warning('The given type is not valid, changing to ''rxns''')
    type = 'rxns';
end

distMat = zeros(n,n);

for i=1:n-1
    row_i = zeros(1,n);
        % Get the adjacency matrix
    adj_i = createRxnConnectivityGraph(models{i}, blackList, biomass);
    % Get the distance matrix that compares the first neighbourhood
    dist_i = sparse(getN1Dist(adj_i));
    parfor j=i+1:n
        adj_j = createRxnConnectivityGraph(models{j}, blackList, biomass);
        % Get the distance matrix that compares the first neighbourhood
        dist_j = sparse(getN1Dist(adj_j));
        % The distance matrices still contain all reactions of both models,
        % whereas here only the intersection should be compared.
        if isequal(type, 'rxns')
            [~, idx_i, idx_j] = intersect(models{i}.rxns, models{j}.rxns);
            % Calculate The R_v coefficient
            row_i(j) = 1 - RvCoefficient(dist_i(idx_i,idx_i), dist_j(idx_j,idx_j) );
        else
            % More pre-processing necessary
            [~, idx_i, idx_j] = intersect(models{i}.EC, models{j}.EC);
            % Calculate The R_v coefficient
            row_i(j) = 1 - RvCoefficient( dist_i(idx_i,idx_i), dist_j(idx_j,idx_j) );
        end
    end
    distMat(i,:) = row_i;
end

% Copy the upper right triangle to the lower left triangle to obtain a
% symmetric distance matrix
for i=1:n
    distMat(:,i) = distMat(i,:);
end

end