function Rv = RvCoefficient(m1, m2)
% Compute the R_v coefficient to determine the congruence of two matrices
% m1 and m2.

if size(m1,1)~=size(m2,1) || size(m1,2)~=size(m2,2)
    error('The number of rows and columns of m1 and m2 have to match.')
end

% If matrices are rectangular, transform them into positive semi-definite matrices
if size(m1,1)~=size(m1,2)
    m1 = m1*m1.';
    m2 = m2*m2.';
end

Rv = trace(m1*m2) / sqrt( trace(m1^2) * trace(m2^2) );

end




