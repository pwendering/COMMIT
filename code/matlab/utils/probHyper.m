function p = probHyper(N, n, K, k)
% calculate the probability for drawing k successes in a hypergeometric 
% distribution defined by a total number of N observations, a test set K 
% (e.g. differentially expressed), number of k successes in the test set K and
% a total number of n successes.
% Input:
%       double N:       total number of observations
%       double n:       total number of successes in N
%       double K:       size of test set
%       double k:       number of successes in the test set


p = (nchoosek(n, k) * nchoosek(N-n, K-k)) / nchoosek(N, K);
end