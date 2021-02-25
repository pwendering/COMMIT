function enrichment = metOntologyEnrichment(metaboliteList,...
    referenceList, ontologyGraph, idList)
%% metOntologyEnrichment: Performs ontology enrichent analysis for metabolites
% Make sure the namespaces for all inputs match!
% Input:
%       cell metaboliteList:        array containing metabolite IDs to
%                                   test for
%       cell referenceList:         reference population of metabolite
%                                   IDs
%       char ontologyGraph:         adjacency matrix of the DAG for
%                                   ontology terms (edges point to the
%                                   root)
%       cell idList:                array containing identifiers of all
%                                   ontology terms in the graph
% Output:
%       table enrichment:           table with columns: "ID", "ontology
%                                   term", "p_value"

warning('off')

if nargin == 4
    if ~iscell(metaboliteList)
        error('metaboliteList must be a cellstr array')
    elseif ~iscell(referenceList)
        error('referenceList must be a cellstr array')
    elseif ~isnumeric(ontologyGraph)
        error('ontologyGraph must be an adjacency matrix')
    elseif size(ontologyGraph,1) ~= size(ontologyGraph,2)
        error('the adjacency matrix is not square')
    elseif ~iscell(idList)
        error('idList must be a cellstr array')
    end
else
    error('Wrong number of input parameters')
end

n = numel(metaboliteList);
N = numel(referenceList);

%% quantify the abundance of each ontology term (115742 for ChEBI (10.01.2020)
abundances_test = getAbundances(ontologyGraph, metaboliteList, idList);
abundances_ref = getAbundances(ontologyGraph, referenceList, idList);

%% Fischer's exact test

% only consider terms that appear in the test set
idx_nz = find(abundances_test);
P = zeros(numel(idx_nz), 1);

% hypergeometric test
for i=1:numel(idx_nz)
    k = abundances_test(idx_nz(i));
    K = abundances_ref(idx_nz(i));
    
    if K<k
        K=k;
    end
    
    % one-sided p-value for upper tail
    P(i) = hygecdf(k, N, K, n, 'upper');
end

enrichment = array2table(P, 'VariableNames', {'p_value'},...
    'RowNames', idList(idx_nz));

end
