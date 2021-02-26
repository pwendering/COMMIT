function identical = compareGPR(model_1, rxn_1, model_2, rxn_2)
%% Compares the GPR rules of two reactions
% (Jaccard index of genes multiplied with Jaccard of logical operators
% Input:
%       struct model_1, model_2:            metabolic models
%       cell/char/double rxn_1, rxn_2:      reactions given as ID or index
% Output:
%       double identical:                   degree of identity
%

%% Check input
if nargin < 4
    error('You did not give enough input arguments')
elseif ~isstruct(model_1) || ~isstruct(model_2)
    error('The model variables have to be of type struct')
end

if ~isnumeric(rxn_1)
    try
        rxn_1 = find(strcmp(model_1.rxns, rxn_1));
    catch
        error('Please give a valid type for the rxn input')
    end
end

if ~isnumeric(rxn_2)
    try
        rxn_2 = find(strcmp(model_2.rxns, rxn_2));
    catch
        error('Please give a valid type for the rxn input')
    end
end

%% Compare gene IDs
% indices of genes
gene_idx_1 = cellfun(@str2double,regexp(model_1.rules{rxn_1}, '\d+', 'match'));
gene_idx_2 = cellfun(@str2double,regexp(model_2.rules{rxn_2}, '\d+', 'match'));

% gene IDs
genes_1 = model_1.genes(gene_idx_1);
genes_2 = model_2.genes(gene_idx_2);

% Output the Jaccard index divided by 2 if operators do not match
rule_1 = model_1.rules{rxn_1};
rule_2 = model_2.rules{rxn_2};

operators_1 = regexp(rule_1, '&|\|', 'match');
operators_2 = regexp(rule_2, '&|\|', 'match');

if numel(operators_1)<1|| numel(operators_2)<1
    operators_factor = 1;
else
    % Jaccard index of operators
    operators_factor = numel(intersect(operators_1, operators_2)) / numel(union(operators_1, operators_2));
end

% Jaccard index of genes muliplied with Jaccard index of operators
identical = numel(intersect(genes_1, genes_2)) / numel(union(genes_1, genes_2)) * operators_factor;

if isnan(identical) || isinf(identical)
    identical = 0;
end
end