% Run conditional gap filling on individual models
tic

disp('----------------------------------------------------------------------')
disp('START')
disp('----------------------------------------------------------------------')

disp('Loading required data...')

% load options script
options;

% load models 
load(modelFile) 

if exist('merged_models', 'var')
    models = merged_models; clear merged_models
end


model_ids = cellfun(@(x)strtok(x.id, '_'), models, 'UniformOutput', false);

% gap-filling database
load(dbFile)

% medium
load(mediumFile)

%{
get the list of complementary media
auxo_media = dir(fullfile(mediaDir, '*.tsv'));
auxo_media = fullfile({auxo_media.folder}, {auxo_media.name})';

create cell arrays for individual auxotrophic media
for j=1:numel(auxo_media)
    tmp_tab = importdata(auxo_media{j});
    tmp_tab = tmp_tab.textdata(2:end,1);
    auxo_media{j} = strcat(translateIDs(tmp_tab, 'met', [], 'ModelSEED',...
        'MNXref', false), '[e]');
end
%}

% workspace containing sequence similarity weights
load(seq_sim_workspace)
idx = ismember(col_labels, model_ids);
seq_sim_mat = seq_sim_mat(:, idx);
col_labels = col_labels(idx);
genes = genes(idx);

n = numel(models);
GF = cell(n,1);
exc = cell(n,1);
bio = zeros(n,1);
gf = zeros(n,1);

disp('-----------------------------------')

fprintf('Starting gap filling of %d models\n\n', n)
for i=1:n
    fprintf('######\t%s\n', model_ids{i})
    
    model = models{i};
%     medium = auxo_media{i};

    % remove possible uptake reactions
    uptake_rxns = (sum(model.S==1,1)==1) & (sum(model.S~=0) == 1);
    model = removeRxns(model, model.rxns(uptake_rxns));
    model = removeRxns(model, model.rxns(contains(model.rxns, 'EX_')));
    % add medium to the model
    model = addExchangeRxn(model, medium);
    % test if model is already functional
    v = cplexlp(-model.c, [], [], model.S, model.b, model.lb, model.ub);
    opt = v(model.c==1);
    
    if opt>=10E-6
        disp('Model already functional')
    end; clear v opt
    
    idx_model = strcmp(model_ids(i), col_labels);
    dbModel_MNXref_balanced.scores = seq_sim_mat(:,idx_model);
    dbModel_MNXref_balanced.genes = repmat({''}, numel(dbModel_MNXref_balanced.rxns),1);
    
    % run fastGapFilling
    [~, ~, GF{i}, addedRxns] = ...
        condFastGapFilling(model, dbModel_MNXref_balanced, [],...
        weights, epsilon, include_sink, true);
    
    % make the stoichiometric matrix sparse
    GF{i}.S = sparse(GF{i}.S);
    options = optimoptions('linprog', 'Display', 'none');
    v = cplexlp(-GF{i}.c, [], [], GF{i}.S, GF{i}.b, GF{i}.lb, GF{i}.ub);
    
    bio(i) = v(logical(GF{i}.c)); clear v
    gf(i) = numel(addedRxns)-sum(contains(addedRxns, 'sink_'));
    
    if include_sink
        % sink reactions were already included in the gap filling
        % LP
        exchanged_mets = regexp(addedRxns, 'MNXM\d+\[.\]', 'match');
        exchanged_mets = [exchanged_mets{:}];
    else
        % LP to find sink reactions that do not decrease the
        % biomass production by a factor
        exchanged_mets = findPotentialExcMets(GF{i}, 0.9);
    end
    
    exchanged_mets = setdiff(exchanged_mets, medium);
    tmp_exc = regexprep(exchanged_mets, 'c', 'e');
    [~, ia] = setdiff(tmp_exc, medium, 'stable');
    
    % get the indices of uptake reactions
    exchange_rxns = findExchangeReactions(GF{i});
    
    mets_added_exc = GF{i}.mets(any(GF{i}.S(:,logical(exchange_rxns)),2));
    mets_added_exc = setdiff(mets_added_exc, [{'BIOMASS[c]'}; medium]);
    
    % set of potentially exchanged metabolites
    P = union(exchanged_mets(ia), mets_added_exc); clear exchanged_mets tmp_exc
    
    
    % extend the set of exchanged metabolites
    exc{i} = P;
    
    disp('---------------')
end

% save the models to workspace
save([modelFile, '_gf'], 'GF', 'exc', 'bio', 'gf')

disp('-----------------------------------')
fprintf('Finished, time: %.2f min.\n', toc)
