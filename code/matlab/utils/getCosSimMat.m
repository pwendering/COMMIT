function similarityMatrix = getCosSimMat(model)
% Compoute the cosine similarity matrix between all columns of a
% stoichiometric matrix
% Input:        struct model:                   metabolic model
% output:       double similarity_matrix:       contains cosine similarity
%                                               values
stoich_mat = full(model.S);
n = numel(model.rxns);
similarity_matrix = zeros(n_rxns);

parfor i=1:n-1
    A = stoich_mat(:,i);
    row_tmp = zeros(1, n);
    for j=i+1:n
        B = stoich_mat(:,j);
        row_tmp(j) = cosineSimilarity(A, B);
    end
    similarity_matrix(i,:) = row_tmp;
end
end