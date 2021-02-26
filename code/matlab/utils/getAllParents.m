function parents = getAllParents(G,V, exclude)
%% parents = getAllParents(G,V, exclude)
% Find all parent nodes in a graph for a given set of nodes.
% The graph is traversed by DFS, returns the indices of the traversed
% nodes.
%  Input:
%           double G:           adjacennnncy matrix of the graph
%           double V:           array containing the vertex indices
%           double exclude:     neighborhood of the vertices to exclude
%                               from the set of all parents
%  Output:
%           cell parents:       array containing the indices of all parent 
%                               nodes for each node respectively (including
%                               the node itself)

if nargin < 3
    exclude = [];
end

parents = cell(numel(V),1);

%% perform DFS for each term and save all traversed nodes
for i=1:numel(V)
    v = V(i);
    parents{i} = graphtraverse(G,v, 'Method', 'DFS');
end

if ismember(exclude, 1:numel(V))
    for i=1:numel(V)
        v = V(i);
        N_exclude = graphtraverse(G,v,'Method', 'BFS', 'Depth', exclude);
        parents{i} = setdiff(parents{i}, N_exclude);
    end
end

end

