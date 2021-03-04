function [GF, EX, gf_order, solutions, exc, gf, bio, dep] = iterativeGapFilling(models, medium, auxo_media, DB, weights, epsilon,...
    include_sink, order, iterations, seq_similarity)
%% iterativeGapFilling
% Function for gap-filling members of a microbial community. Each member is
% used as a seed in 100 random successsions each iteratively
% gap-filling the models in the community. Whenever metabolites in the
% preceding models are exchanged, these metabolites are added to
% the gap filling medium (exchange reactions in the model) for the
% following community members. The final solution is the minimal
% solution over all iterations and all seeds.
% Input:
%       cell models:            array containing the metabolic models as
%                               struct objects
%       cell medium:            array containing the metabolites of the
%                               minimal medium (make sure the same namespace is
%                               used)
%       cell auxo_media:        array containing predicted auxotrophic media for each model
%                               individually; if empty, the minimal medium
%                               will be used
%       struct DB:              universal database that is used for gap
%                               filling (contains fields 'rxns', 'mets', 'S', 'transport', 'lb', 'ub'
%       struct weights:         contains weights for the different
%                               types of reactions that can be added
%                               from the database: 'transport',
%                               'metabolic', 'model' (reactions in
%                               original model) and 'reverse'
%                               (reactions that have the opposite
%                               direction to irreversible reactions
%                               in the original model)
%       double epsilon:         threshold for the biomass reaction
%       logical include_sink:   if true, the gap filling objective will
%                               include potentially exported metabolites as
%                               sink reactions
%       double order:           abundances of the organisms or OTUs
%                               represented by the individual models given
%                               as absulute values or percentages
%       double iterations:      number of iterations
%       char seq_similarity:    name of the workspace where sequence
%                               similarity weights and genes are stored
%                               'seq_sim_mat, 'genes', 'col_labels'
% Output:
%       cell GF:                array containing the gap-filled metabolic
%                               models
%       cell EX:                array containing the exchanged metabolites
%       double gf_order:        arrax containing the optimal gap filling
%                               order used for generatin GF
%       double solutions:       matrix containing all random successions
%       double exc:             matrix containing the numbers of exchangable
%                               metabolites per model model per iteration
%       double gf:              matrix containing the numbers of added
%                               reactions per model per iteration
%       double bio:             matrix containing all biomass fluxes for
%                               each model in every iteration
%       double dep:             fractions of compounds in the first auxo
%                               medium that can be produced by subsequent models
warning('off', 'all') % show no warnings
format shortg % for time format

%% Parameters and Initialization
if nargin < 1
    error('USAGE: iterativeGapFilling(models, medium, DB)')
elseif nargin < 3
    error('You did not give enough input arguments; USAGE: iterativeGapFilling(models, medium, DB')
elseif ~isstruct(models{1}) && ~iscellstr(models)
    error('Please check the type of your ''models'' input variable')
elseif ~iscellstr(medium)
    error('Please check the type of your ''medium'' input variable')
elseif ~isstruct(DB) || ~all(isfield(DB, {'rxns', 'mets', 'S', 'transport', 'lb', 'ub'}))
    error('Please check if the database is of type struct and if it contains all required fields')
elseif isempty(auxo_media)
    auxo_media = repmat({medium}, numel(models), 1);
elseif exist('order', 'var') && ~isempty(order)
    if numel(order)==numel(models)
        [~,gf_order] = sort(order);
        iterations = 0;
    else
        error('Your ''order'' argument has incorrect length')
    end
end

if exist('seq_similarity', 'var') && ischar(seq_similarity) && ~isempty(seq_similarity)
    fprintf('\nReading sequence similarity file...\n');
    load(seq_similarity, 'seq_sim_mat', 'genes', 'col_labels')
    
    % take subset of models in the current community
    if iscellstr(models)
        model_ids = cellfun(@(x)regexp(x, '\w+$', 'match'), models);
    else
        model_ids = cellfun(@(x)strtok(x.id, '_'), models,...
            'UniformOutput', false);
    end
    
    idx = ismember(col_labels, model_ids);
    genes = genes(idx);
    seq_sim_mat = seq_sim_mat(:, idx);
    col_labels = col_labels(idx);
    
    seq_sim.matrix = seq_sim_mat;
    seq_sim.genes = genes;
    seq_sim.labels = col_labels;
    clear seq_sim_mat genes col_labels
else
    seq_sim = [];
    DB.scores = ones(numel(DB.rxns), 1);
    DB.genes = repmat({''}, numel(DB.rxns), 1);
end

% initial time point
t_init = clock;
% numer of models to be gap filled
n = numel(models);
% default solution (optimal succession for gap filling)
solutions = zeros(iterations, n);
% array to save the number of added reactions per sequence
gf = zeros(iterations,n);
% array to save the number of exchanged metabolites
exc = zeros(iterations,n);
% array to save the sums of biomasses of the individual successions
bio = zeros(iterations,n);
% array to save the fraction of metabolites in the auxotrophic medium of
% the first model can be exported by subsequent models
dep = zeros(iterations, 1);

%% Iterate over all models (seeds)
disp('----------------------------------------------------------------------')
disp('START')
disp('----------------------------------------------------------------------')

if iterations > 1
    fprintf('Testing %d random successions for %d models\n\n', iterations, n)
end

parfor i=1:iterations
    warning('off', 'all') % show no warnings
    changeCobraSolver('matlab','LP');

    format shortg % for time format
    t_1 = clock;
    
    % create random successsion and add it to candidate set
    solutions(i,:) = randperm(n);
    
    [gf(i,:), exc(i,:), bio(i,:), dep(i)] = gapFillSuccession(models, solutions(i,:), DB, medium, auxo_media, seq_sim,...
        weights, epsilon, include_sink, false);
    
    t_2 = clock;
    t_diff = t_2 - t_1;
    t_diff = abs(t_diff(4))*60+t_diff(5)+t_diff(6)/60; % minutes
    
    fprintf(['Iteration #%d\n'...
        'Time:\t%.2f min\n',...
        'Total number of added reactions:\t\t%.f\n',...
        'Total number of exchanged metabolites:\t\t%.f\n',...
        'Sum of optimal biomass fluxes:\t\t\t%.2f\n\n'],...
        i, t_diff, sum(gf(i,:)), sum(exc(i,:)), sum(bio(i,:)))
end

%% repeat iterative gap filling for the optimal solution
if iterations~=0
    
    % calculate row sums
    g = sum(gf, 2);
    b = sum(bio, 2);
    e = sum(exc, 2);
    
    disp('----------------------------------------------------------------------')
    
    % first find the iteration(s) with minimal number of added reactions
    idx_opt = find(g==min(g));
    
    % if solution not unique, take the solution with the highest prpportion
    % of produced seed medium components
    if numel(idx_opt) > 1
        idx_opt = idx_opt(dep(idx_opt)==max(dep(idx_opt)));
    end
    
    % if solution not unique, take the solution with the highest numeber of
    % exchanged metabolites amoung the remaining solutions
    if numel(idx_opt) > 1
        idx_opt = idx_opt(e(idx_opt)==max(e(idx_opt)));
    end
    
    % if solution not unique, take the solution with the highest sum of
    % biomasses amoung the remaining solutions
    if numel(idx_opt) > 1
        idx_opt = idx_opt(b(idx_opt)==max(b(idx_opt)));
    end
    
    gf_order = solutions(idx_opt, :);
    
    if sum(idx_opt) > 1
        % if still multiple results, choose one at random
        gf_order = gf_order(randperm(size(gf_order,1),1), :);
    end
    
end

disp('Iterative Gap filling using the optimal succession')

% array containing exchanged metabolites
EX = {};
% gap filling solution
GF = cell(n,1);

M = auxo_media{gf_order(1)};

for i = gf_order
    
    % if models are given as a cell array, select the current model or
    % load is from workspace if file names are given
    if iscellstr(models)
        load(models{i})
    else
        model = models{i};
    end
    
    % remove possible uptake reactions
    uptake_rxns = (sum(model.S==1,1)==1) & (sum(model.S~=0)==1);
    model = removeRxns(model, model.rxns(uptake_rxns));
    model = removeRxns(model, model.rxns(contains(model.rxns, 'EX_')));
    
    % add the medium as exchange reactions
    model = addExchangeRxn(model, M);
    
    if ~isempty(seq_sim)
        % add respective sequence similarity values and genes
        tmp_idx = find(strcmp(strtok(model.id, '_'), seq_sim.labels));
        DB.scores = seq_sim.matrix(:, tmp_idx);
        DB.genes = seq_sim.genes{tmp_idx};
        clear tmp_idx
    end
    
    % conditional FastGapFilling
    [~, ~, gfModel, addedRxns] = ...
        condFastGapFilling(model, DB, EX, weights, epsilon, include_sink, false);
    clear model
    
%     v = cplexlp(-gfModel.c, [], [], gfModel.S, gfModel.b, gfModel.lb, gfModel.ub);
    v = optimizeCbModel(gfModel);
    v = v.x;
    
    if include_sink
        % sink reactions were already included in the gap filling
        % LP
        exchanged_mets = regexp(addedRxns, 'MNXM\d+\[.\]', 'match');
        exchanged_mets = [exchanged_mets{:}];
    else
        % LP to find sink reactions that do not decrease the
        % biomass production by a factor
        exchanged_mets = findPotentialExcMets(gfModel, 0.9);
    end
    
    % indices of metabolites unique to exchanged metabolites
    ext_exc = regexprep(exchanged_mets, '\[c\]', '[e]');
    ext_exc = setdiff(ext_exc, M, 'stable');
    
    % get all external metabolites that take part in active reactions but
    % are not part of the medium
    ext_mets = gfModel.mets(any(gfModel.S(:, v > epsilon),2)); clear v
    ext_mets = setdiff(ext_mets(contains(ext_mets, '[e]')), M);
    
    % set of potentially exchanged metabolites
    P = vertcat(ext_exc, ext_mets); clear ext_mets
    
    % add transport reactions for every cytosolic metabolite that
    % can be exchanged
    cyt_exc = exchanged_mets(contains(exchanged_mets, '[c]'));
    gfModel = addExportReactions(gfModel, cyt_exc);
    
    % sink reactions for all exchangeable metabolites
    [~,gfModel] = evalc(['addSinkReactions(gfModel, P,',...
        'zeros(size(P)), repmat(1000, size(P)));']);
    
    % extend the set of exported metabolites
    EX = union(EX, P);
    
    GF{i} = gfModel;
    
    % minimal medium will be used for all subsequent gap fillings
    M = medium;
end

EX = EX';
t_2 = clock;
t_diff = t_2 - t_init;
t_diff = abs(t_diff(4))*60+t_diff(5)+t_diff(6)/60; % minutes
t_diff = [round(t_diff/60), mod(t_diff,60)];
fprintf('\nTime:\t%dh %.2fmin\n', t_diff(1), t_diff(2))

end