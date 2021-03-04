function [reaction_sets, solutions, consistModel, addedRxns] = ...
    condFastGapFilling(model, dbModel, excMets, weights, epsilon, include_sink, verbose)
% Runs the FastGapfilling algorithm (Latendresse 2014, BMC) on the
% universal database generated using prepareFastGapFilling. The original
% model and the database model are converted to irreversible by splitting
% all reversible reactions into two irreversible reactions. The reactions
% that are shared between the model and the database are removed from the 
% database model by ID. The reactions from the Non-recongized
% reactions are added via their metabolites keeping the reacion labels of
% the original model. The gap filling algorithm itself is an LP that
% maximizes a a trade-off between optimal biomass production and insertion
% of additional reactions. This LP is solved in a binary search with
% decreasing weight on the biomass reaction. The consistent matrix is then
% extracted and forms the new stoichiometric matrix of the metabolic model.
% 
% Input:
%           struct model:               metabolic model to be gap-filled
%           struct dbModel:             stoichiometric model of universal
%                                       database obtained from prepareFastGapFilling
%           cell excMets:               array containing metabolites whose
%                                       exchange is allowed to be added,
%                                       i.e. special weights can be
%                                       assigned in the weights struct 
%           struct weights:             contains weights for the different
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
%           logical include_sink:       if true, the gap filling objective will
%                                       include potentially exported metabolites as
%                                       sink reactions
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
if ~exist('verbose', 'var')
    verbose = true;
end

tic
%% Prepare the datbase matrix 
% convert to database model to irreversible
dbModel_irr = convertModelToIrreversible(dbModel);
% genes field
rev = intersect(find(dbModel.lb), find(dbModel.ub));
if isfield(dbModel, 'genes')
     dbModel_irr.genes = vertcat(dbModel_irr.genes, dbModel_irr.genes(rev));
end
clear dbModel

% convert the input model to irreversible
model_irr = convertModelToIrreversible(model);
rxns_model_irr = model_irr.rxns;
biomass_id = model.rxns(logical(model.c));
biomass_idx = find(logical(model.c));
biomass_sink_id = {'sink_BIOMASS[c]'};
biomass_sink_idx = find(strcmp(model.rxns, biomass_sink_id));

% Get the ids of all irreversible reactions, they can be assigned different 
% weights to because they can be reversed (indices will be the same as in
% the original model because reverse reactions are just appended)
irrev = union(find(model.lb==0), find(model.ub==0));
rev_cand_idx = setdiff(irrev, [biomass_idx, biomass_sink_idx]);
rev_cand = setdiff(rxns_model_irr(irrev), [biomass_id, biomass_sink_id]);
clear model irrev

if verbose
    fprintf('\nPreparing the universal database...\n')
end

% general set difference of reactions in the database and the model
idx_matched = ismember(strtok(dbModel_irr.rxns, '_r'),...
    strtok(rxns_model_irr, '_r'));
dbModel_irr.rxns(idx_matched) = [];
dbModel_irr.lb(idx_matched) = [];
dbModel_irr.ub(idx_matched) = [];
dbModel_irr.transport(idx_matched) = [];
dbModel_irr.S(:,idx_matched) = [];
if isfield(dbModel_irr, 'genes')
    dbModel_irr.genes(idx_matched) = [];
    dbModel_irr.scores(idx_matched) = [];
end
if isfield(dbModel_irr, 'rxnPermeability')
    dbModel_irr.rxnPermeability(idx_matched) = [];
end
clear idx_matched

% add metabolites that are not present in the database
mets_not_in_db = model_irr.mets(~ismember(model_irr.mets, dbModel_irr.mets));
% metabolite IDs
dbModel_irr.mets = vertcat(dbModel_irr.mets,...
    mets_not_in_db);
% metabolite names
if isfield(dbModel_irr, 'metNames')
    dbModel_irr.metNames = vertcat(dbModel_irr.metNames,...
        mets_not_in_db);
end
% stoichiometric matrix
dbModel_irr.S = [dbModel_irr.S; zeros(numel(mets_not_in_db),...
    size(dbModel_irr.S,2))];
% permeability vector
if isfield(dbModel_irr, 'metPermeability')
    dbModel_irr.metPermeability = vertcat(dbModel_irr.metPermeability,...
    zeros(numel(mets_not_in_db), 1));
end
clear mets_not_in_db

% get the indices of exchange reactions in the database model
exchange_db = findExchangeReactions(dbModel_irr);

% index translation for metabolites in the model to the database
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

if include_sink
    % add sink reactions for all metabolites
    S_sink = ones(size(dbModel_irr.S,1), 1);
    % stoichiometric coefficients
    S_sink = -1*S_sink;
    % exclude metabolites that are not likely to diffuse
    S_sink(~logical(dbModel_irr.metPermeability), :) = 0;
    
    % find transport reactions, add metabolites that have a transport reaction
    % in the model between cytosol and extracellular space
    transported = {};
    for i=1:numel(model_irr.rxns)
        
        tmp_mets = model_irr.mets(any(model_irr.S(:,i),2));
        tmp_comps = regexp(tmp_mets, '\[.+?\]', 'match');
        tmp_comps = unique([tmp_comps{:}]);
        
        if numel(tmp_comps)>1
            % if transport reaction, find the transported metabolite(s)
            tmp_mets_ID = strtok(tmp_mets, '[');
            tmp_occurrence = sum(string(tmp_mets_ID)==string(tmp_mets_ID'));
            transported = vertcat(transported, tmp_mets_ID(tmp_occurrence>1));
            clear tmp_mets tmp_mets_ID tmp_occurrence tmp_comps
        end
    end
    % find the indices of the transported metabolites in the database model
    transported = strcat(unique(transported), '[c]');
    tmp_met_idx = cellfun(@(x)find(strcmp(dbModel_irr.mets,x)),...
        transported, 'UniformOutput', false);
    tmp_met_idx = [tmp_met_idx{:}];
    % add sink reaction for each transported metabolite
    
    S_sink(tmp_met_idx) = -1;
    % indices of metabolites with sink reactions (from both permeability
    % and transport reactions)
    tmp_met_idx = any(S_sink,2);
    % create the matrix and delete empty columns
    S_sink = sparse(diag(S_sink));
    S_sink = S_sink(:, any(S_sink,1));
    % metabolite names for the sink reaction IDs
    tmp_mets = dbModel_irr.mets(tmp_met_idx);
else
    S_sink = [];
    tmp_mets = {};
end; clear tmp_met_idx

% create a matrix of exchange reactions for the provided list of allowed
% exchange reactions
if ~isempty(excMets) && isfield(weights, 'uptake')
    
    % empty matrix
    S_ex = sparse(size(dbModel_irr.S,1), numel(excMets));
    
    % in case that the exchanged metabolites contain compartments
    % identifiers, change all to 'e'
    excMets = strtok(excMets, '[');
    excMets = strcat(excMets, '[e]');
    
    % find indices for all metabolites in the list
    idx_exc = cellfun(@(x)find(strcmp(x, dbModel_irr.mets)), excMets,...
        'UniformOutput', false);
    
    excMets = excMets(~cellfun('isempty', idx_exc));
    
    for i=1:numel(excMets)
        S_ex(idx_exc{i}, i) = 1;
    end
    
    
elseif ~isempty(excMets) && ~isfield(weights, 'uptake')
        warning('No weight for ''uptake'' found in the provided weights')
else
    S_ex = [];
end



% add model and reversible candidates to the database
n_rxns_db = size(dbModel_irr.S,2);
% stoichiometric matrix
dbModel_irr.S = [dbModel_irr.S, sparse(S_add), sparse(S_rev_cand), S_sink, S_ex];
% reactions
dbModel_irr.rxns = vertcat(dbModel_irr.rxns, rxns_model_irr,...
    strcat(rxns_model_irr(rev_cand_idx), '_r'),...
    strcat('sink_',tmp_mets),...
    strcat('EX_', excMets));
% lower boundary
dbModel_irr.lb = vertcat(dbModel_irr.lb, model_irr.lb,...
    model_irr.lb(rev_cand_idx), zeros(size(S_sink,2), 1),...
    zeros(size(S_ex,2), 1));
% upper boundary
dbModel_irr.ub = vertcat(dbModel_irr.ub, model_irr.ub,...
    model_irr.ub(rev_cand_idx), repmat(1000, size(S_sink,2),1),...
    repmat(1000, size(S_ex,2),1));
% transport
dbModel_irr.transport = vertcat(dbModel_irr.transport,...
    zeros(numel(rxns_model_irr) + numel(rev_cand_idx) + size(S_sink,2)...
    + size(S_ex,2), 1));
% genes
if isfield(dbModel_irr, 'genes')
    dbModel_irr.genes = vertcat(dbModel_irr.genes,...
        repmat({''}, numel(rxns_model_irr) + numel(rev_cand_idx) + size(S_sink,2)...
        + size(S_ex,2), 1));
    % e-values / scores
    dbModel_irr.scores = vertcat(dbModel_irr.scores,...
        ones(numel(rxns_model_irr) + numel(rev_cand_idx) + size(S_sink,2)...
        + size(S_ex,2), 1));
end
if isfield(dbModel_irr, 'rxnPermeability')
    dbModel_irr.rxnPermeability = vertcat(dbModel_irr.rxnPermeability,...
        zeros(numel(rxns_model_irr) + numel(rev_cand_idx) + size(S_sink,2)...
        + size(S_ex,2),1));
end

% indices for assignment of weights
rxns_model = n_rxns_db+1:n_rxns_db+numel(rxns_model_irr);
rev_cand_idx = max(rxns_model)+1:max(rxns_model) + numel(rev_cand_idx);
rxns_sink = max(rev_cand_idx)+1:max(rev_cand_idx) + size(S_sink,2);
rxns_ex = size(dbModel_irr.S,2)-size(S_ex, 2)+1:size(dbModel_irr.S,2);

clear S_add S_rev_cand S_sink n_rxns_db idx_translation tmp_mets transported tmp_met_idx

tmp_toc = toc;
if verbose
    fprintf('\n\t> finished adding reactions (%.0fs)\n', toc)
end


%% Run the FastGapFilling LP with binary search

% define b vector for steady state
beq = zeros(size(dbModel_irr.S, 1), 1);

% index of the biomass reaction
biomass = strcmp(dbModel_irr.rxns, biomass_id);

% remove exchange reactions from the gap filling matrix
dbModel_irr.S(:, logical(exchange_db)) = 0;

% define sequence similarity weights
t = 10E-6;
if isfield(dbModel_irr, 'scores')
    seq_evidence = dbModel_irr.scores < t;
else
    seq_evidence = [];
end

% determine the number of reactions in the database and the model
m = size(dbModel_irr.S, 2);

% define the objective for the LPs
f = zeros(m,1);
f(logical(dbModel_irr.transport)) = weights.transport;
f(~logical(dbModel_irr.transport)) = weights.metabolic;
f(rev_cand_idx) = weights.reverse;
f(seq_evidence) = weights.sequence;
if isfield(dbModel_irr, 'rxnPermeability')
    f(logical(dbModel_irr.rxnPermeability)) = weights.permeable;
end
f(rxns_model) = weights.model;
f(rxns_sink) = weights.sink;
f(logical(exchange_db)) = weights.exchange;
f(rxns_ex) = weights.uptake;

% initial reaction set
reaction_sets = {};

% initial solutions set for the flux vectors
solutions = {};

% initial values for binary search
alpha = 0;
beta = 2*m;
beta = 10000;
% lower and upper bounds
lb = dbModel_irr.lb;
ub = dbModel_irr.ub;
ub(biomass) = 2.81;
% Start binary search
precision = 10E-6;
if verbose
    fprintf('\nStarting the binary search for the gap filling LP...\n')
end

while abs(alpha - beta) > 1
    % Weighting factor for biomass reaction
   delta = floor(mean([alpha, beta]));
   
   % re-define the objective for the biomass
   f(biomass) = -delta;
   
   % Solve the LP
   solution = cplexlp(f, [], [], dbModel_irr.S, beq, lb, ub);
%    options = optimoptions('linprog', 'Display', 'none');
%    solution = linprog(f, [], [], dbModel_irr.S, beq, lb, ub, options);

   if solution(biomass) >= epsilon
       % consider the reactions that have not been in the model before
       nz = solution >= precision;
       % v vector of reaction fluxes
       solutions{end+1} = solution(union(find(nz),rxns_model));
       nz(rxns_model) = 0;
       % reaction that would be added by this solution
       reaction_sets{end+1} = dbModel_irr.rxns(nz);
       
       
       if numel(reaction_sets)==1 || numel(reaction_sets{end}) <= numel(reaction_sets{end-1})
           beta = delta;
       else
           alpha = delta;
       end
       
       % will be used as the final solution
       S = solution;
   else
       alpha = delta;
   end

end

if verbose
    fprintf('\n\t> finished binary search (%.0fs)\n', toc-tmp_toc)
end
tmp_toc = toc;

clear Aeq nz

%% Update the input model
if ~isempty(reaction_sets)
    if verbose
        fprintf('\nGap filling successful, adding reactions to model...\n')
    end
    
    % start with the irreversible model (so fields irrespecitive of gap
    % filling are kept)
    consistModel = model_irr;
    
    % find the active reactions in the gap filling solution
    nz = S >= precision;
    
    % flux through the biomass reaction
    flux_biomass = S(biomass);
    
    % take the last (smallest) set of added reactions
    addedRxns = reaction_sets{end};
    
    % number of added sink reactions
    n_sink = sum(contains(addedRxns, 'sink'));
    
    % extract the stoichiometric matrix from the gap-filling solution (all
    % added reactions, no sink reactions and all reactions that were alreacy
    % contained in the model
    nz(rxns_sink) = 0;
    nz(rxns_model) = 1;
    consistModel.S = full(dbModel_irr.S(:, nz));
    
    % take only the rows that have entries
    met_idx_add = any(consistModel.S, 2);
    consistModel.S = consistModel.S(met_idx_add, :);
    
    % reactions of the consistent model
    consistModel.rxns = dbModel_irr.rxns(nz);
    
    % total number of added reactions
    n_added_rxns = numel(addedRxns);
    
    % number of reversed reactions
    n_reversed = contains(addedRxns, rev_cand);
    
    % number of transport reactions containing metabolites with high
    % permeability
    n_high_permeability = sum(ismember(dbModel_irr.rxns(logical(dbModel_irr.rxnPermeability)), addedRxns));
    
    % number of added transport reactions
    nz(rxns_model) = 0;
    n_added_trans = sum(dbModel_irr.transport(nz));
    
    % number of added exchange reactions
    n_added_exchange = sum(nz(logical(exchange_db)));
    
    % number of added allowed exchange reactions
    n_added_uptake = sum(nz(rxns_ex));
    
    % number of reactions with sequence support
    if isfield(dbModel_irr, 'scores')
        n_sequence_support = sum(dbModel_irr.scores(nz)<t);
    else
        n_sequence_support = 0;
    end
    
    
    % update all reaction-related fields:
    idx_new_rxns = cell2mat(cellfun(@(x)find(strcmp(x, consistModel.rxns)), addedRxns,...
        'UniformOutput', false));
    
    % boundaries
    nz(rxns_model) = 1;
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
 
    % E.C. numbers
    if isfield(model_irr, 'EC')
        consistModel.EC = cellfun(@(x)model_irr.EC(strcmp(x, rxns_model_irr)),...
            consistModel.rxns, 'UniformOutput', false);
        consistModel.EC(idx_new_rxns) = {''};
    end

    % add new genes and according rules to the model
    if isfield(dbModel_irr, 'genes')
        nz(rxns_model) = 0;
        new_genes = dbModel_irr.genes(nz); % same order as reactions

        
        for i=1:numel(new_genes)
            
            if ~ismember(new_genes(i), consistModel.genes) && ~isempty(new_genes{i})
                % if the gene is not contained, add the gene and the
                % new index for the GPR rule is the number of genes
                consistModel.genes = vertcat(consistModel.genes, new_genes(i));
                if isfield(consistModel, 'geneNames')
                    consistModel.geneNames = vertcat(consistModel.geneNames, new_genes(i));
                end
                new_idx = num2str(numel(consistModel.genes));
                % if the gene is contained, find the index
            elseif isempty(new_genes{i})
                new_idx = 0;
            else
                new_idx = num2str(find(strcmp(consistModel.genes, new_genes{i})));
            end
            
            if new_idx~=0
                % only simple one-gene rules are added
                new_rule = {strcat('x(', new_idx, ')')};
                % add the new rule to the merged model
                consistModel.rules = vertcat(consistModel.rules, new_rule);
            else
                % if no gene associated to the reaction, add an empty entry
                consistModel.rules = vertcat(consistModel.rules, {''});
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
    idx_added_mets_db = cell2mat(cellfun(@(x)find(strcmp(x, dbModel_irr.mets)),...
        addedMets, 'UniformOutput', false));
    consistModel.metNames(idx_new_mets) = dbModel_irr.metNames(idx_added_mets_db);
    
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
        consistModel.metCharges(idx_tmp) = model_irr.metCharges;
        clear idx_tmp
    end
    
    % metabolite formulae
    if isfield(model_irr, 'metFormulas')
        consistModel.metFormulas = repmat({''}, numel(model_irr.mets), 1);
        idx_tmp = cell2mat(cellfun(@(x)find(ismember(consistModel.mets, x)),...
            model_irr.mets, 'UniformOutput', false));
        consistModel.metFormulas(idx_tmp) = model_irr.metFormulas;
        clear idx_tmp
    end
    
    if isfield(consistModel, 'grRules')
        consistModel = rmfield(consistModel, 'grRules');
    end
    
    if verbose
        fprintf('\n\t> finished adding reactions to the model (%.0fs)\n\n', toc-tmp_toc)
        fprintf('\n\t> %d reaction(s) in total have been added\n', n_added_rxns)
        fprintf('\t> %d of the added reactions have sequence support\n', n_sequence_support)
        fprintf('\t> %d transport reaction(s) have been added\n', n_added_trans)
        fprintf('\t\t> %d transport reaction(s) include a metabolite with high permeability\n', n_high_permeability)
        fprintf('\t> %d exchange reaction(s) have been added\n', n_added_exchange + n_added_uptake)
        fprintf('\t> %d of them are allowed uptake reaction(s)\n', n_added_uptake)
        if include_sink
            fprintf('\t> %d sink reactions have been included for highly-permeable metabolites\n', n_sink)
        end
        fprintf('\t> %d metabolic reaction(s)\n', n_added_rxns-n_added_trans-n_sink-n_added_uptake)
        fprintf('\t> %d reaction(s) have been made reversible\n', sum(n_reversed))
        fprintf('\t> %d metabolite(s) have been added to the model.\n', numel(addedMets))
        fprintf('\n\tFlux through biomass reaction:\t %.2f\n', flux_biomass)
    end
else
    warning('Gap filling did not result in a feasible solution')
    consistModel = struct;
    addedRxns = {};
end
if verbose
    fprintf('\nTotal time: %.0fs\n', toc); 
end
end
