function ec = ecNumbersFromModelSEED(reactionIDs, ModelSEEDRxnEC)
% Obtain the Enzyme Commission (EC) numbers that are assigned to reactions in the ModelSEED database
% Input:
%           cell reactionIDs:           contains the ModelSEED reaction ids
%           table ModelSEEDReactionDB:  table with two columns
%                                       (1: reaction IDs; 2: E.C. numbers)
% Output:
%           cell ec:                    contains the EC numbers

if ~iscell(reactionIDs)
    error('The reaction IDs are not given in a cell array.')
end
% Read the file ModelSEEDReactionDB

n = numel(reactionIDs);
ec = cell(n,1);

for i=1:n
    rxn = regexp(reactionIDs{i}, 'rxn[0-9]*', 'match');
    if ~isempty(rxn)
        idx = strcmp(ModelSEEDRxnEC{:,1}, rxn);
        if any(idx)
            if ~any(strcmp(ModelSEEDRxnEC{idx, 2}, 'null'))
                ec(i) = ModelSEEDRxnEC{idx, 2};
            else
                ec(i) = {''};
            end
        end
    else 
        ec(i) = {''};
    end
end
end



