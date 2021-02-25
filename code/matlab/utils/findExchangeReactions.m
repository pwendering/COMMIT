function ex = findExchangeReactions(model)
% Find exchange reactions in a metabolic model
%
sink_rxns = (sum(model.S==-1,1)==1) & (sum(model.S~=0) == 1);
uptake_rxns = (sum(model.S==1,1)==1) & (sum(model.S~=0) == 1);
ex = sink_rxns + uptake_rxns;
end