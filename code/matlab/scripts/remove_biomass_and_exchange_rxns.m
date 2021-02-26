% remove exchange reactions and biomass reactions
options
clearvars -except topDir
habitats = {'Soil', 'Leaf', 'Root'};
modelTopDir = fullfile(topDir, 'data', 'models');

% KBase models
disp('KBase')
modelDir = fullfile(modelTopDir, 'kbase');
for i=1:numel(habitats)
    disp(['>', habitat])
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_genes_translated.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        ec = model.EC;
        idx_bio = strcmp(model.rxns, 'BIOMASS_Reaction');
        
        % find exchange reactions
        exchange_rxns = logical(findExchangeReactions(model))';
        
        to_remove = model.rxns(logical(exchange_rxns+idx_bio));
        model = removeRxns(model, to_remove, 'metFlag', false);
        
        % deal with EC numbers separately since there can be errors if the
        % number of metabolites and reactions is equal
        model.EC = ec(~logical(exchange_rxns + idx_bio));
        
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    save(workspace, 'models')
    clear models
end

% RAVEN models
disp('RAVEN')

modelDir = fullfile(modelTopDir, 'raven',' HMMer10E-50');
for i=1:numel(habitats)
    disp(['>', habitat])
    workspace = fullfile(modelDir, strcat(habitats{i}, '_converted.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        ec = model.EC;
        idx_bio = strcmp(model.rxns, 'BIOMASS_Reaction');
        
        % find exchange reactions
        exchange_rxns = logical(findExchangeReactions(model))';
        
        to_remove = model.rxns(logical(exchange_rxns+idx_bio));
        model = removeRxns(model, to_remove, 'metFlag', false);
        
        % deal with EC numbers separately since there can be errors if the
        % number of metabolites and reactions is equal
        model.EC = ec(~logical(exchange_rxns + idx_bio));
        
        models{j} = model;
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    save(workspace, 'models')
    clear models
end


% CarveMe
disp('CarveMe')
modelDir = fullfile(modelTopDir, 'carveme');
for i=1:numel(habitats)
    disp(['>', habitat])
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        ec = model.EC;
        idx_bio = strcmp(model.rxns, 'BIOMASS_Reaction');
        
        % find exchange reactions
        exchange_rxns = logical(findExchangeReactions(model))';
        
        to_remove = model.rxns(logical(exchange_rxns+idx_bio));
        model = removeRxns(model, to_remove, 'metFlag', false);
        
        % deal with EC numbers separately since there can be errors if the
        % number of metabolites and reactions is equal
        model.EC = ec(~logical(exchange_rxns + idx_bio));
        
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    save(workspace, 'models')
    clear models
    
end

% AuReMe
disp('AuReMe')
modelDir = fullfile(modelTopDir, 'aureme');
for i=1:numel(habitats)
    disp(['>', habitat])
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_genes_translated.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        ec = model.EC;
        idx_bio = strcmp(model.rxns, 'BIOMASS_Reaction');
        
        % find exchange reactions
        exchange_rxns = logical(findExchangeReactions(model))';
        
        to_remove = model.rxns(logical(exchange_rxns+idx_bio));
        model = removeRxns(model, to_remove, 'metFlag', false);
        
        % deal with EC numbers separately since there can be errors if the
        % number of metabolites and reactions is equal
        model.EC = ec(~logical(exchange_rxns + idx_bio));
        
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    save(workspace, 'models')
    clear models
end
clear topDir