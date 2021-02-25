function ec_pruned = pruneEC(ec, level)
% Prune a list of EC numbers to a specified level.
% Input:
%           cell ec:            Array containing the ec numbers
%           double level:       Number that specifies the level on
%                               which the EC numbers should be compared
% Output:
%           cell ec_pruned:     Array that contains the pruned EC numbers

if ~iscellstr(ec)
    error('Input must be of type cellstr')
end

ec_pruned = cellfun(@(x)strsplit(x, '.'),ec,'UniformOutput', false);
ec_pruned = cellfun(@(x)strjoin(x(1:level),'.'), ec_pruned, 'UniformOutput', false);
end