function excIdx = findExchangeReactions(model)
%% exchangeIdx = findExchangeReactions(model)
% Find exchange reactions in a metabolic model
% Input:
%       struct model:           metabolic model with field 'S'
% Output:
%       logical excIdx:         boolean array indicating indices of
%                               exchange reactions
sinkRxnIdx= (sum(model.S==-1,1)==1) & (sum(model.S~=0) == 1);
uptakeRxnIdx = (sum(model.S==1,1)==1) & (sum(model.S~=0) == 1);
excIdx = sinkRxnIdx + uptakeRxnIdx;
end