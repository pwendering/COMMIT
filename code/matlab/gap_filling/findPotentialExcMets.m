function transportable = findPotentialExcMets(model, alpha)
% Find transport and sink reactions for permeable metabolites that can 
% be added to the model without decreasing the biomass production 
% to a large extend 
% Input:
%       struct model:           metabolic model
% Output:
%       cell transportable:     metabolites for which transport and
%                               exchange reactions likely
%% Parameters
epsilon = 1E-5;
if nargin < 2
    alpha = 0.9;
end

%% Get permeability for all metabolites in the model
permeable_idx = logical(getPermeabilityWeights(model.mets));
mets_permeable = strtok(model.mets(permeable_idx));

% only cytosolic metabolites
mets_permeable = mets_permeable(contains(mets_permeable, '[c]'));

%% find transported metabolites
transported = {};
for i=1:numel(model.rxns)
    tmp_mets = model.mets(any(model.S(:,i),2));
    tmp_comps = regexp(tmp_mets, '\[.+?\]', 'match');
    tmp_comps = unique([tmp_comps{:}]);
    if numel(tmp_comps)>1
        % transport reaction
        tmp_mets_ID = strtok(tmp_mets, '[');
        tmp_occurrence = sum(string(tmp_mets_ID)==string(tmp_mets_ID'));
        transported = vertcat(transported, tmp_mets_ID(tmp_occurrence>1));
        clear tmp_mets tmp_mets_ID tmp_occurrence tmp_comps
    end
end

% find the indices of the transported and permeable metabolites in the model
transported = unique([strcat(unique(transported), '[c]'); mets_permeable]);
tmp_met_idx = cellfun(@(x)find(strcmp(model.mets,x)),...
    transported, 'UniformOutput', false);
tmp_met_idx = [tmp_met_idx{:}]';
tmp_mets = model.mets(tmp_met_idx);

%% Find the optimal biomass value
model = convertModelToIrreversible(model);
% solution = optimizeCbModel(model);
% solution = solution.x;
solution = cplexlp(-model.c, [], [], model.S, model.b, model.lb, model.ub);
opt = solution(logical(model.c));

%% Sink reactions
S_sink = sparse(eye(size(model.S,1)));
% add sink reaction for each transported metabolite
S_sink = -1*S_sink;
% remove empty columns
S_sink = S_sink(:, tmp_met_idx);

%% Linear program 
% solve a linear program to include as many transport reaction for
% highly-permeable metabolites as possible without decreasing the biomass
% by a factor alpha

% combine stoichiometric matrix and sink reactions
Aeq = [model.S, S_sink];
% b vector for steady state
beq = zeros(size(Aeq,1), 1);

% objective (maximize the flux through the sink reactions)
f = [zeros(size(model.S,2),1); zeros(size(S_sink,2),1)];
f(size(model.S,2)+1:end) = -1;

% define lower and upper bounds for the LP
lb = [model.lb; zeros(size(S_sink,2),1)];
lb(find(model.c)) = alpha*opt;
ub = [model.ub; repmat(1000, size(S_sink,2),1)];

%{
% solve the LP
lp.A = Aeq;
lp.b = beq;
lp.lb = lb;
lp.ub = ub;
lp.c = -f;
lp.osense = -1;
lp.csense = repmat('E',1,size(lp.A,1));
solution = solveCobraLP(lp);
solution = solution.full;
%}

solution = cplexlp(f, [], [], Aeq, beq, lb, ub);

%% Find transportable metabolites
% active sink reaction above epsilon
active_sink = solution(size(model.S,2)+1:end)>=epsilon;

transportable = tmp_mets(active_sink);

end