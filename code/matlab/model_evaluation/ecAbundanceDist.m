function distMat = ecAbundanceDist(models, allECNumbers)
% Creates distance matrix based on the abundance of every Enzyme Commission
% (EC) number present in the model. The resulting vectors are compared using the
% Spearman rank correlation, which is substracted from 1 to obtain a
% distance betwenn 0 and 2.
% Input:
%               cell models:        array that contains the metabolic
%                                   models, all model structs should at least have the 'ec' field
%               cell allECNumbers:  array that contains a list of all EC
%                                   numbers
%
% Output:   double distMat:         distance matrix

n = numel(models);
m = numel(allECNumbers);

for i=1:n
    if ~isfield(models{i}, 'EC')
        error('At least one model does not have a field ec containing Enzyme Commission numbers')
    end
end

if ~iscellstr(allECNumbers)
    error('Please give the list of all EC numbers as a cellstr')
end

% Regex pattern for EC numbers that contain three or more levels
pattern = '[0-9A-Z]+\.[0-9A-Z]+\.[0-9A-Z]+\.[0-9A-Z]+';

% Calculate the abundance for every EC number in every model
abundance = cell(n, 1);
parfor i=1:n
    % Obtain a list of all ec numbers in the model that contain all four
    % levels since the list of all EC numbers only contains complete
    % numbers
    ec = cellfun(@(x)regexp(x,pattern,'match'),...
        models{i}.EC,'UniformOutput', false);
    ec = [ec{:}];
    adj = 1 - numel(models{i}.rxns) / m;
    
    % Calculate the abundance of exact matches
    res = cellfun(@(x)strcmp(x, ec), allECNumbers, 'UniformOutput', false);
    abundance{i} = cellfun(@(x)sum(x), res)*adj;
end
clear models

% Calculate the distances
distMat = zeros(n,n);

for i=1:n-1
    for j=i+1:n
        distMat(i,j) = 1 - corr(abundance{i}, abundance{j},...
            'type', 'Spearman');
        distMat(j,i) = distMat(i,j);
    end
end

end