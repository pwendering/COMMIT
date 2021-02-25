function [reaction_sets, solutions, consistModel, addedRxns] = FastGapFilling(model, dbModel, weights, epsilon)
% Runs the FastGapfilling algorithm (Latendresse 2014, BMC) on the
% universal database generated using prepareFastGapFilling. The original
% model and the database model are converted to irreversible by splitting
% all reversible reactions into two irreversible reactions. The reactions
% that are shared between the model and the database are checked for
% identity so the labels match the correct direction. Non-recongized
% reactions are added via their metabolites keeping the reacion labels of
% the original model. The gap filling algorithm itself is an LP that
% maximizes a a trade-off between optimal biomass production and insertion
% of additional reactions. This LP is solved in a binary search with
% decreasing weight on the biomass reaction. The consistent matrix is then
% extracted and forms the new stoichiometric matrix of the metabolic model.
% 
% Input:
%           struct model:               metabolic model to be gap-filled
%           struct databaseModel:       stoichiometric model of universal
%                                       database obtained from prepareFastGapFilling
%           struct weights:              contains weights for the different
%                                       types of reactions that can be added
%                                       from the database: 'transport',
%                                       'metabolic', 'model' (reactions in
%                                       original model), 'reverse'
%                                       (reactions that have the opposite
%                                       direction to irreversible reactions
%                                       in the original model), 'exchange'
%                                       (for exchange reactions) and
%                                       'sequence' (vector: weights based on
%                                       sequence support)
%           double epsilon:             threshold for the biomass reaction
% Output:
%           cell reaction_sets:         all reaction sets with f_biomass > epsilon
%                                       that result from the binary search
%                                       (as cellstr)
%           cell solutions:             all solution vectors with f_biomass > epsilon
%                                       (as double)
%           struct consistModel:        model with updated stoichiometric
%                                       matrix including the smallest (or last) 
%                                       gap filling reaction set (all
%                                       reacitons irreversible)
%           cellstr addedRxns:          array that contains all reactions
%                                       that were added to consistModel
tic
%% Prepare the datbase matrix 
% convert to database model to irreversible
dbModel_irr = convertModelToIrreversible(dbModel);

% genes field
rev = intersect(find(dbModel.lb), find(dbModel.ub));
dbModel_irr.genes = vertcat(dbModel_irr.genes, dbModel_irr.genes(rev));
clear dbModel

% convert the input model to irreversible
model_irr = convertModelToIrreversible(model);
rxns_model_irr = model_irr.rxns;
biomass_id = model.rxns(logical(model.c));
biomass_idx = find(logical(model.c));
biomass_sink_id = {'EX_BIOMASS_c'};
biomass_sink_idx = find(strcmp(model.rxns, 'EX_BIOMASS_c'));
% Get the ids of all irreversible reactions, they can be assigned different 
% weights to because they can be reversed (indices will be the same as in
% the original model because reverse reactions are just appended)
irrev = union(find(model.lb==0), find(model.ub==0));
rev_cand_idx = setdiff(irrev, [biomass_idx, biomass_sink_idx]);
rev_cand = setdiff(rxns_model_irr(irrev), [biomass_id, biomass_sink_id]);
clear model irrev
% model_irr.ub(model_irr.ub>0) = 2000;
% model_irr.lb(model_irr.lb<0) = -2000;

fprintf('\nPreparing the universal database...\n')

% general set difference of reactions in the database and the model
idx_matched = ismember(strtok(dbModel_irr.rxns, '_r'),...
    strtok(rxns_model_irr, '_r'));
dbModel_irr.rxns(idx_matched) = [];
dbModel_irr.lb(idx_matched) = [];
dbModel_irr.ub(idx_matched) = [];
dbModel_irr.transport(idx_matched) = [];
dbModel_irr.S(:,idx_matched) = [];
dbModel_irr.genes(idx_matched) = [];
dbModel_irr.scores(idx_matched) = [];
clear idx_matched

% add metabolites that are not present in the database
mets_not_in_db = model_irr.mets(~ismember(model_irr.mets, dbModel_irr.mets));
% metabolite IDs
dbModel_irr.mets = vertcat(dbModel_irr.mets,...
    mets_not_in_db);
% stoichiometric matrix
dbModel_irr.S = [dbModel_irr.S; zeros(numel(mets_not_in_db),...
    size(dbModel_irr.S,2))];
% permeability vector
dbModel_irr.permeable_idx = vertcat(dbModel_irr.rxnPermeability,...
    zeros(numel(mets_not_in_db), 1));
clear mets_not_in_db
% build index translation for metabolites in the model to the database
idx_translation = cell2mat(cellfun(@(x)find(strcmp(x, dbModel_irr.mets)),...
    model_irr.mets, 'UniformOutput', false));

S_add = zeros(size(dbModel_irr.S, 1), size(model_irr.S, 2));

% expand the stoichiometric matrix of the model to the matrix of the
% database
for i=1:numel(rxns_model_irr)
    % get the metabolite indices for the current reaction in the model
    met_idx = find(model_irr.S(:,i));
    % metabolite indices in the database
    met_idx_db = idx_translation(met_idx);
    % stoichiometric coefficients
    coeff = model_irr.S(met_idx,i);
    % add coefficients to the the reaction in S_add
    S_add(met_idx_db, i) = coeff;
end

% add reactions with opposite directions to originally irreversible reactions in the model 
S_rev_cand = zeros(size(dbModel_irr.S, 1), numel(rev_cand));
for i=1:numel(rev_cand_idx)
    rxn_idx = rev_cand_idx(i);
    % get the metabolite indices for the current reaction in the model
    met_idx = find(model_irr.S(:,rxn_idx));
    % metabolite indices in the database
    met_idx_db = idx_translation(met_idx);
    % stoichiometric coefficients
    coeff = model_irr.S(met_idx,rxn_idx);
    % add coefficients to the the reaction in S_add
    S_rev_cand(met_idx_db, i) = -coeff;
end

% add model and reversible candidates to the database
n_rxns_db = size(dbModel_irr.S,2);
% stoichiometric matrix
dbModel_irr.S = [dbModel_irr.S, sparse(S_add), sparse(S_rev_cand)];
% reactions
dbModel_irr.rxns = vertcat(dbModel_irr.rxns, rxns_model_irr,...
    strcat(rxns_model_irr(rev_cand_idx), '_r'));
% lower boundary
dbModel_irr.lb = vertcat(dbModel_irr.lb, model_irr.lb,...
    model_irr.lb(rev_cand_idx));
% upper boundary
dbModel_irr.ub = vertcat(dbModel_irr.ub, model_irr.ub,...
    model_irr.ub(rev_cand_idx));
% transport
dbModel_irr.transport = vertcat(dbModel_irr.transport,...
    zeros(numel(rxns_model_irr)+numel(rev_cand_idx),1));
% genes
dbModel_irr.genes = vertcat(dbModel_irr.genes,...
    repmat({''}, numel(rxns_model_irr)+numel(rev_cand_idx), 1));
% e-values / scores
dbModel_irr.scores = vertcat(dbModel_irr.scores,...
    ones(numel(rxns_model_irr)+numel(rev_cand_idx), 1));

% indices for assignment of weights
rxns_model = n_rxns_db+1:n_rxns_db+numel(rxns_model_irr);
rxns_rev_cand = numel(rxns_model_irr)+1:size(dbModel_irr.S,2);
clear S_add S_rev_cand n_rxns_db idx_translation

tmp_toc = toc;
fprintf('\n\t> finished adding reactions (%.0fs)\n', toc)


%% Run the FastGapFilling LP with binary search
% generate the stoichiometric matrix for database and model
Aeq = dbModel_irr.S;

% define b vector for steady state
beq = zeros(size(dbModel_irr.S, 1), 1);

% index of the biomass reaction
biomass = strcmp(dbModel_irr.rxns, biomass_id);

% get the indices of exchange reactions in the database model
sink_rxns = (sum(dbModel_irr.S==-1,1)==1) & (sum(dbModel_irr.S~=0) == 1);
uptake_rxns = (sum(dbModel_irr.S==1,1)==1) & (sum(dbModel_irr.S~=0) == 1);
exchange_rxns = sink_rxns + uptake_rxns;

% define sequence similarity weights
t = 10E-6;
seq_evidence = dbModel_irr.scores < t;

% Upper and lower boundaries
lb = dbModel_irr.lb;
ub = dbModel_irr.ub;

% determine the number of reactions in the database and the model
m = size(dbModel_irr.S, 2);
n = size(model_irr.S, 2);

% define the objective for the LPs
f = zeros(m,1);
f(logical(dbModel_irr.transport)) = weights.transport;
f(~logical(dbModel_irr.transport)) = weights.metabolic;
f(logical(exchange_rxns)) = weights.exchange;
f(rev_cand_idx) = weights.reverse;
f(logical(dbModel_irr.permeable_idx)) = weights.permeable;
f(seq_evidence) = weights.sequence;
f(rxns_model) = weights.model;
% initial reaction set
reaction_sets = {};

% initial solutions set for the flux vectors
solutions = {};

% initial values for binary search
alpha = 0;
beta = n + m;

% Start binary search
precision = 10E-6;
fprintf('\nStarting the binary search for the gap filling LP...\n')
while abs(alpha - beta) > 1
    % Weighting factor for biomass reaction
   delta = floor(mean([alpha, beta]));
   
   % re-define the objective for the biomass
   f(biomass) = -delta;
   
   % Solve the LP
   solution = cplexlp(f, [], [], Aeq, beq, lb, ub);
   
   if solution(biomass) >= epsilon
       % consider the reactions that have not been in the model before
       nz = solution > precision;
       % v vector of reaction fluxes
       solutions{end+1} = solution(union(find(nz),rxns_model));
       nz(rxns_model) = 0;
       % reaction that would be added by this solution
       reaction_sets{end+1} = dbModel_irr.rxns(nz);
       % will be used as the final solution
       S = solution;

       beta = delta;
   else
       alpha = delta;
   end

end
fprintf('\n\t> finished binary search (%.0fs)\n', toc-tmp_toc)
tmp_toc = toc;

clear Aeq nz

%% Update the input model
if ~isempty(reaction_sets)
    fprintf('\nGap filling successful, adding reactions to model...\n')
    consistModel = model_irr;
    % Subsect the database matrix with the reactions from the original model
    % and the added reactions
    nz = S > precision;
    nz(rxns_model) = 1;
    consistModel.S = full(dbModel_irr.S(:, nz));
    met_idx_add = any(consistModel.S, 2);
    consistModel.S = consistModel.S(met_idx_add, :);
    consistModel.rxns = dbModel_irr.rxns(nz);
    addedRxns = reaction_sets{end};
    reversed = contains(addedRxns, rev_cand);
    n_high_permeability = sum(ismember(dbModel_irr.rxns(logical(dbModel_irr.permeable_idx)), addedRxns));
    n_added_rxns = numel(addedRxns);
    nz(rxns_model) = 0;
    n_added_trans = sum(dbModel_irr.transport(nz));
    n_sequence_support = sum(dbModel_irr.scores(nz)<t);
    nz(rxns_model) = 1;
    
    
    % update all reaction-related fields:
    idx_new_rxns = cell2mat(cellfun(@(x)find(strcmp(x, consistModel.rxns)), addedRxns,...
        'UniformOutput', false));
    
    % boundaries
    consistModel.lb = dbModel_irr.lb(nz);
    consistModel.ub = dbModel_irr.ub(nz);
    
    % reaction names
    consistModel.rxnNames = cellfun(@(x)model_irr.rxnNames(strcmp(x, rxns_model_irr)),...
        consistModel.rxns, 'UniformOutput', false);
    consistModel.rxnNames(idx_new_rxns) = {''};
    
    % reaction notes
    consistModel.rxnNotes = cellfun(@(x)model_irr.rxnNotes(strcmp(x, rxns_model_irr)),...
        consistModel.rxns, 'UniformOutput', false);
    consistModel.rxnNotes(idx_new_rxns) = {'gf'};
    
    % objective ('c')
    consistModel.c = zeros(size(consistModel.S, 2), 1);
    consistModel.c(strcmp(consistModel.rxns, biomass_id)) = 1;
    
    % subsystems
    consistModel.subSystems = cellfun(@(x)model_irr.subSystems(strcmp(x, rxns_model_irr)),...
        consistModel.rxns, 'UniformOutput', false);
    consistModel.subSystems(idx_new_rxns) = {''};
    
    % gene rules
    consistModel.rules = cellfun(@(x)model_irr.rules(strcmp(x, rxns_model_irr)),...
        consistModel.rxns, 'UniformOutput', false);
    consistModel.rules(idx_new_rxns) = {''};
    
    % E.C. numbers
    consistModel.EC = cellfun(@(x)model_irr.EC(strcmp(x, rxns_model_irr)),...
        consistModel.rxns, 'UniformOutput', false);
    consistModel.EC(idx_new_rxns) = {''};
    
    % add new genes and according rules to the model
    if isfield(model_irr, 'genes')
        nz(rxns_model) = 0;
        new_genes = dbModel_irr.genes(nz);
        n_genes = numel(consistModel.genes);
        consistModel.genes = vertcat(consistModel.genes, new_genes);
        consistModel.geneNames = vertcat(consistModel.geneNames, new_genes);
        for i=1:numel(new_genes)
            if ~isempty(new_genes{i})
                model.rules{idx_new_rxns(i)} = ...
                    strcat('x(', num2str(n_genes+i), ')');
            end
            
        end
    end
    
    % get the metabolites associated with these reactions
    consistModel.mets = dbModel_irr.mets(met_idx_add);

    % update all metabolite-related fields:
    addedMets = setdiff(consistModel.mets, model_irr.mets);
    idx_new_mets = cell2mat(cellfun(@(x)find(strcmp(x, consistModel.mets)), addedMets,...
        'UniformOutput', false));
    
    % metabolite names
    consistModel.metNames = cellfun(@(x)model_irr.metNames(strcmp(x, model_irr.mets)),...
        consistModel.mets, 'UniformOutput', false);
     consistModel.metNames(idx_new_mets) = {''};
    
    % b vector
    consistModel.b = zeros(size(consistModel.S, 1), 1);
    
    % constraint type ('csense')
    consistModel.csense = char(cellfun(@(x)model_irr.csense(strcmp(x, model_irr.mets)),...
        consistModel.mets, 'UniformOutput', false));
    consistModel.csense(idx_new_mets) = 'E';
    
    % metabolite charge
    if isfield(model_irr, 'metCharges')
        consistModel.metCharges = nan(numel(consistModel.mets), 1);
        consistModel.metCharges = cell2mat(cellfun(@(x)model_irr.metCharges(strcmp(x, model_irr.mets)),...
            consistModel.mets, 'UniformOutput', false));
        idx_tmp = cell2mat(cellfun(@(x)find(ismember(consistModel.mets, x)),...
            model_irr.mets, 'UniformOutput', false));
        consistModel.metCharges(idx_tmp) = cell2mat(cellfun(@(x)model_irr.metCharges(strcmp(x, model_irr.mets)),...
            consistModel.mets, 'UniformOutput', false));
        clear idx_tmp
    end

    fprintf('\n\t> finished adding reactions to the model (%.0fs)\n\n', toc-tmp_toc)
    
    fprintf('\n\t> %d reaction(s) in total have been added\n', n_added_rxns)
    fprintf('\t> %d of the added reactions have sequence support\n', n_sequence_support)
    fprintf('\t> %d transport or exchange reaction(s)\n', n_added_trans)
    fprintf('\t> %d transport reaction include a metabolite with high permeability\n', n_high_permeability)
    fprintf('\t> %d metabolic reaction(s)\n', n_added_rxns-n_added_trans)
    fprintf('\t> %d reaction(s) have been made reversible\n', sum(reversed))
    fprintf('\t> %d metabolite(s) have been added to the model.\n', numel(addedMets))
    
else
    warning('Gap filling did not result in a feasible solution')
    consistModel = struct;
    addedRxns = {};
end

fprintf('\nTotal time: %.0fs\n', toc); 
end
