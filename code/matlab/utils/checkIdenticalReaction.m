function [identical, reversible] = checkIdenticalReaction(model_1, rxn_1, model_2, rxn_2, blackList, check_opposite)
% Check for the identity of two reactions in two different models. Does not
% account for false assignment of stoichiometric coefficients to the
% metabolites
% Input:
%           struct model_1, model_2:        metabolic models that contain
%                                           the fields 'S', 'mets' and 'rxns'
%           char or double rxn_1, rxn_2:    reaction ids that should be
%                                           compared (or given as inidices)
%           cellstr blackList:              list of metabolites that should
%                                           be excluded from the comparison 
%                                           (e.g. H+ is oftentimes not
%                                           considered)
%           logical checkOpposite:          whether irreversible reactions
%                                           that share the same metabolites 
%                                           and coefficients should be
%                                           reported as identical
% Output:
%           logical identical:              1 if reactions are identical, 0
%                                           otherwise
%           logical reversible:             1 if reactions have opposite
%                                           directions (and are identical);
%                                           0 otherwise
if nargin < 5
    blackList = {};
    check_opposite = 0;
end
if nargin < 6
    check_opposite = 0;
end

% Threshold for identity
epsilon = 1E-6;
t = 1 - epsilon;

%% find the coefficients for the first reaction:
% row and column indices and exclude metabolites from blackList
if ~isnumeric(rxn_1)
    idx_rxn_1 = strcmp(model_1.rxns, rxn_1);
else
    idx_rxn_1 = rxn_1;
end

nz_rxn_1 = any(model_1.S(:,idx_rxn_1), 2);

% find indices of metabolites to exclude and set them to 0
nz_rxn_1(ismember(model_1.mets, blackList)) = 0;

% get metabolite ids and coefficients
coeff_1 = model_1.S(nz_rxn_1, idx_rxn_1);
mets_1 = model_1.mets(nz_rxn_1);
% sort the coefficients with respect to alphabetic order of metabolites
[mets_1, idx_1] = sort(mets_1);
coeff_1 = coeff_1(idx_1);

%% find the coefficients for the second reaction:
% row and column indices and exclude metabolites from blackList
if  ~isnumeric(rxn_2)
    idx_rxn_2 = strcmp(model_2.rxns, rxn_2);
else
    idx_rxn_2 = rxn_2;
end
nz_rxn_2 = any(model_2.S(:, idx_rxn_2), 2);


% find indices of metabolites to exclude and set them to 0
nz_rxn_2(ismember(model_2.mets, blackList)) = 0;

% get metabolite ids and coefficients
coeff_2 = model_2.S(nz_rxn_2, idx_rxn_2);
mets_2 = model_2.mets(nz_rxn_2);
% sort the coefficients with respect to alphabetic order of metabolites
[mets_2, idx_2] = sort(mets_2);
coeff_2 = coeff_2(idx_2);

%% test for identity
identical = 0;
reversible = 0;
if isequal(mets_1, mets_2)
    similarity = cosineSimilarity(coeff_1, coeff_2);
    if ~check_opposite
        identical = similarity >= t;
    elseif similarity <= -t
        identical = 1;
        reversible = 1;
    elseif similarity >= t
        identical = 1;
    end
end

identical = logical(identical);
end