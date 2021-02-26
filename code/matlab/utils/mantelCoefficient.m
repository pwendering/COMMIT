function r_M = mantelCoefficient(A,B)
%% r_M = mantelCoefficient(A,B)
% Calculates the Mantel coefficient for two n x n distance matrices A and B

if nargin == 2
    if ~isequal(size(A), size(B)) || ~isequal(size(A,1), size(A,2))
        error('The matrices A and B must have the same dimensions')
    elseif ~isnumeric(A) || ~isnumeric(B)
        error('Both matrices have to be numeric')
    end
else
    error('Wrong number of input arguments')
end

% set the diagonal values to 0
A = A - diag(A)'*eye(size(A));
B = B - diag(B)'*eye(size(B));

% Row and column dimension
n = size(A,2);

% mean values of the off-diagonal values
a = sum(sum(A)) / ( n * n - n );
b = sum(sum(B)) / ( n * n - n );

% substract the mean from all values and take only the upper right triangle
A_diff = triu(A - repmat(a, n), 1);
B_diff = triu(B - repmat(b, n), 1);

% covariance of upper off-diagonal entries 
cov = trace(A_diff * B_diff');

% pooled standard deviation of A and B upper off-diagonal entries
pool_sd = sqrt(...
    trace(A_diff * A_diff') * ...
    trace(B_diff * B_diff')...
    );

r_M = cov / pool_sd;

end