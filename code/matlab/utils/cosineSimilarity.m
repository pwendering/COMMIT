function similarity = cosineSimilarity(A, B)
%% similarity = cosineSimilarity(A, B)
% Calculate the cosine similarity by 
% similarity = ( A .* B ) / ( ||A|| * ||B|| )
% Input:
%       double A, B:        vectors to be compared
% Output:
%       double similarity:  cosine similarity between A and B

if isvector(A) && isvector(B)
    similarity = dot( A, B ) / ( norm(A) * norm(B) );
else
    error('At least one of the inputs is not a vector')
end
end
