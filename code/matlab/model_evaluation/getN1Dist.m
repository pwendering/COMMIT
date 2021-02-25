function distMat = getN1Dist(adjMat)
% Create a distance matrix between all reactions or Enzyme Commission (EC)
% numbers that are connected via metabolites in the rxnAdjMat (obtained from 
% createRxnConnectivityGraph(model, blackList)). The intersection of the 
% first neighbourhoods of the vertices is returned as the 'distance'.
% Input:
%           double rxnAdjMat:       adjacency matrix of the reactions or EC
%                                   numbers
% Output:
%           double distMat:         distance matrix

if size(adjMat,1)~=size(adjMat,2)
    error('The given matrix is not square.')
end

n = size(adjMat,2);
distMat = zeros(n,n);

for i=1:n-1
    for j=i+1:n
        distMat(i,j) = sum(adjMat(:,i).*adjMat(:,j));
        distMat(j,i) = distMat(i,j);
    end
end 

end