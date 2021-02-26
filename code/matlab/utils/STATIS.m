function [S, q] = STATIS(M)
%% [S, q] = STATIS(M)
% Computes a matrix consensus using the STATIC method (Lavit et al. 1994, 
% Computational Statistics & Data Analysis)
% Input:
%       cell M:         array containing the matrices to be merged
% Output:
%       double S:       consensus matrix
%       double q:       ratio of the first eigenvalue and the sum of all
%                       eigenvalues

n = numel(M);

%% Compare the matrices using the Rv-coefficient
C = zeros(n,n);

for i=1:n
    for j=1:n
        C(i,j) = RvCoefficient(M{i}, M{j});
    end
end

%% Find the leading eigenvector and re-scale it so its sum of elements is one
% eigenvectors and eigenvalues of C
[V,D] = eig(C);
% sort the eigenvalues in descending order
[d, ia] = sort(diag(D), 'descend'); clear D
% re-arrange the eigenvector matrix
V = V(:,ia); clear ia
% select the eigenvector corresponding to the largest eigenvalue
a = V(:,1); clear V
% scale the selected eigenvector
a = a / sum(a);
% divide the first eigenvalue by the sum of eigenvalues
q = d(1) / sum(d);

%% compute the compromise matrix as the weighted sum of all matrices in M
S = zeros(size(M{1}));

for i=1:n
    S = S + a(i) * M{i};
end

end
