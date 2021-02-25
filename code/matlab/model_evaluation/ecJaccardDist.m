function distMat = ecJaccardDist(models, level)
% Creates a distance matrix based on the Jaccard distance between two set
% of Enzyme Commission (EC) numbers on a level from 1 to 4.
% Input:
%               cell models:        cell array that contains the metabolic
%                                   models, all model structs should at least have the 'EC' field
%               double level:       Number that specifies the level on
%                                   which the EC numbers should be compared
% Output:   double distMat:         distance matrix

n = numel(models);

for i=1:n
    if ~isfield(models{i}, 'EC')
        error('At least one model does not have a field EC containing Enzyme Commission numbers')
    end
end


% Regex pattern for EC numbers that contain three or more levels
pattern = '[0-9A-Z]+\.[0-9A-Z\-]+\.[0-9A-Z\-\.]+';

distMat = zeros(n,n);

% write 'EC' fields into cell array
EC = cell(n, 1);
for i=1:n
    EC{i} = models{i}.EC;
end
clear models

parfor i=1:n-1
    % Get the pruned EC numbers for model i:
    ec_i = cellfun(@(x)regexp(x,pattern,'match'),...
        EC{i},'UniformOutput', false);
    % Sometimes, EC numbers are given with less levels without replacing with
    % '-', so this character has to be added
    ec_i = cellfun(@(x)completeEC(x),[ec_i{:}], 'UniformOutput', false);
    ec_i = pruneEC(ec_i, level);
    %ec_i =[ec_i{:}];
    
    row_i = zeros(1,n);
    
    for j=i+1:n
        % Get the pruned EC numbers for model j
        ec_j = cellfun(@(x)regexp(x,pattern,'match'),...
            EC{j},'UniformOutput', false);
        ec_j = cellfun(@(x)completeEC(x),[ec_j{:}], 'UniformOutput', false);
        ec_j = pruneEC(ec_j, level);
        %ec_j =[ec_j{:}];
        % Calculate the Jaccard distance
        row_i(j) = 1-numel(intersect(ec_i, ec_j))/numel(union(ec_i, ec_j));
    end
    
    distMat(i,:) = row_i;
end

% Copy the upper right triangle to the lower left triangle to obtain a
% symmetric distance matrix
for i=1:n
    distMat(:,i) = distMat(i,:);
end

end





