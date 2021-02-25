function abundances = getAbundances(G, metList, idList, exclude)
%% getAbundances quantify the abundance of every ontology term in a list
% First, all parents of each term are determined, the abundances of every
% of the resulting terms is quantified
% Input:
%       double G:           adjacency matrix of the ontology DAG (edges point
%                           to the root)
%       cell metList:       array containing metabolite IDs (terms) to query
%                           for
%       cell idList:        list of all ontology terms
%       double exclude:     neighborhood of the vertices to exclude
%                           from the set of all parents
% Output:
%       double abundances:  abundance for every ontology terms respectively

%% find indices of the test set metabolites  ==> vertex indices
vertex_indices = cell2mat(cellfun(@(x)find(strcmp(idList,x)), metList,...
    'UniformOutput', false));

%% get all parents for the terms of interest
ontology_terms = getAllParents(G, vertex_indices, exclude);
ontology_terms = [ontology_terms{:}];

%% translate vertex indices back to ontology terms
ontology_terms = idList(ontology_terms);

%% quantify the abundance of each ontology term
abundances = cell2mat(cellfun(@(x)sum(ismember(ontology_terms, x)),...
    idList, 'UniformOutput', false));

end
