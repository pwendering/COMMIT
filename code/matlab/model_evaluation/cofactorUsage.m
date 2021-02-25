function cofactorMat = cofactorUsage(models, cofactorList)
% Creates a matrix of size #models x #cofactors that displays the number of
% occurrence of the given cofactors, respectively.
% Input:
%               cell models:            array that contains the metabolic models
%               cell cofactorList:      array that contains the cofactors
%                                       the are in the same namespace as 
%                                       the metabolites in the model
% Output:
%               double cofactorMat:     matrix that contains the number
%                                       of occurrences for each cofactor

n = numel(models);
m = numel(cofactorList);

cofactorMat = zeros(n,m);

if iscellstr(cofactorList)
    for i=1:m
        cofactorList(i) = {cofactorList(i)};
    end
end

for i=1:n
    mets = models{i}.mets;
    mets = regexprep(mets, '\[[a-z]\]$', '');
    for j=1:m
        for k=1:numel(cofactorList{j})
            cofactorMat(i,j) = 0;
            % get the row index if the cofactor if present
            idx = strcmp(cofactorList{j}{k},mets);
            if sum(idx)~=0
                % calculate the sum of the row entries that have been converted
                % to logical not to consider stoichiometric coefficients
                cofactorMat(i,j) = cofactorMat(i,j) + length(find(logical(models{i}.S(idx,:))==1));
            end
        end
        
    end
end


end