function adj_mat = createRxnConnectivityGraph(model, blackList, biomass)
% Creates an adjacency matrix of reactions of a metabolic model.
% Two reactions are connected if they share at least one metabolite. 
% Highly abundant metabolites such as e.g H2O or reduction equivalents are 
% removed (blackList) before constructing the graph.
% Input:
%           struct model:           model that contains at least the fields
%                                   'mets' and 'S'
%           cell blacklist:         metabolite ids that should be excluded
%           char biomass:           reaction id of the biomass reaction, if
%                                   no biomass reaction is defined, input 'none'
% Output:   double adj_mat:         adjacency matrix of the resulting graph


if ~isfield(model, 'S') || ~isfield(model, 'mets')
    error('The model does not contain one of the fields S or mets or neither')
end

matrix = model.S;

% remove higly abundant metabolites from the network
remove_idx = cellfun(@(x)ismember(x, blackList), regexprep(model.mets, '\[.\]$', '') ,...
    'UniformOutput', false);
remove_idx = ~logical(cell2mat(remove_idx));
% remove the biomass reaction
if sum(contains(model.rxns, biomass))~=1 && ~isequal(biomass, 'none')
    error('The biomass reaction id is not uniqe')
else
    remove_idx(contains(model.rxns, biomass)) = 0;
end

matrix = matrix(remove_idx,:);

n = size(matrix, 2);
adj_mat = zeros(n,n);

for i=1:n-1
    for j=i+1:n
        adj_mat(i,j) = sum(matrix(:,i).*matrix(:,j)) ~= 0;
        adj_mat(j,i) = adj_mat(i,j);
    end
end

end


