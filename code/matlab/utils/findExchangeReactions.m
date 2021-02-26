function exchangeIdx = findExchangeReactions(model)
% Find exchange reactions in a metabolic model
sinkRxnIdx= (sum(model.S==-1,1)==1) & (sum(model.S~=0) == 1);
uptakeRxnIdx = (sum(model.S==1,1)==1) & (sum(model.S~=0) == 1);
exchangeIdx = sinkRxnIdx + uptakeRxnIdx;
end