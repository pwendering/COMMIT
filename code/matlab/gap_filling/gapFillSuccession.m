function [gf, exc, bio, dep] = gapFillSuccession(models, S, DB, medium, auxo_media, seq_sim,...
    weights, epsilon, include_sink, verbose)
% run conditinoal FastGapFilling with custom weights on a given succession
% Inputs see iterativeGapFilling

% number of models
n = numel(models);

% number of added reactions per model during gap filling
gf = zeros(1,n);

% biomass fluxes of the gap-filled models
bio = zeros(1,n);

% number of potentially exchanged metabolites per gap-filled model in
% succession
exc = zeros(1,n);

% empty set of exchanged metabolites
EX = {};

% gap filling medium
auxo_first = auxo_media{S(1)};
M = auxo_first; clear auxo_media


% fill gaps iteratively

for i = S
    
    % if models are given as a cell array, select the current model or
    % load is from workspace if file names are given
    if iscellstr(models)
        load(models{i})
    else
        model = models{i};
    end
    
    % remove possible uptake reactions
    uptake_rxns = (sum(model.S==1,1)==1) & (sum(model.S~=0) == 1);
    model = removeRxns(model, model.rxns(uptake_rxns));
    model = removeRxns(model, model.rxns(contains(model.rxns, 'EX_')));
    
    % add medium components to the model as exchange reactions if
    % not already present
    model = addExchangeRxn(model, unique(M));
    
    if ~isempty(seq_sim)
        % add respective sequence similarity values and genes
        tmp_idx = find(strcmp(strtok(model.id, '_'), seq_sim.labels));
        DB.scores = seq_sim.matrix(:, tmp_idx);
        DB.genes = seq_sim.genes{tmp_idx};
        clear tmp_idx
    end

    % run conditional FastGapFilling
    [~, ~, gfModel, addedRxns] = ...
        condFastGapFilling(model, DB, EX, weights, epsilon, include_sink, verbose);
    clear model

    % add the optimal biomass value of the current gap filled model
    v = optimizeCbModel(gfModel).x;
    bio(i) = v(logical(gfModel.c));
    
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
    exchanged_mets = regexprep(exchanged_mets, '\[c\]', '[e]');
    exchanged_mets = setdiff(exchanged_mets, M, 'stable');
    
    % get all external metabolites that take part in active extracellular reactions
    % are not part of the medium
    ext_mets = gfModel.mets(any(gfModel.S(:, v > epsilon),2)); clear v
    ext_mets = setdiff(ext_mets(contains(ext_mets, '[e]')), M);
    
    % set of potentially exchanged metabolites
    P = vertcat(exchanged_mets, ext_mets); clear exchanged_mets ext_mets
    
    % extend the set of exchanged metabolites
    EX = union(EX, P);
    
    % increase the counter of added reactions
    gf(i) = numel(addedRxns)-sum(contains(addedRxns, 'sink_'));
    
    exc(i) = numel(P);
    
    if ~iscellstr(models)
        models{i} = 0;
    end
    
    % minimal medium will be used for all subsequent gap fillings
    M = medium;
    
end

% fraction of exported metabolites and the auxotrophic medium for the
% first model (dependency of first model on the others)
% if dep is close to 1, it is very likely, the the first model is largely
% independent of the auxo medium its the required nutrients can be produced
% by other members
dep = numel(intersect(auxo_first, EX)) / numel(auxo_first);

end
