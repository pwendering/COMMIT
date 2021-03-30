% add formulae to models

formulaTab = readtable('data/tables/MNXref/MNXref_MET_FORMULAE.csv',...
    'ReadVariableNames', true);
habitats = {'Soil', 'Leaf', 'Root'};

% KBase models
modelDir = 'data/models/kbase';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    load(workspace)
    for j=1:numel(models)
        models{j} = addMetFormulae(models{j}, formulaTab);
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    save(workspace, 'models')
    clear models
end

% RAVEN models
modelDir = 'data/models/raven/HMMer10E-50';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    load(workspace)
    for j=1:numel(models)
        models{j} = addMetFormulae(models{j}, formulaTab);
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_models_metFormulas.mat'));
    save(workspace, 'models')
    clear models
end

% CarveMe
modelDir = 'data/models/carveme';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    load(workspace)
    for j=1:numel(models)
        models{j} = addMetFormulae(models{j}, formulaTab);
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    save(workspace, 'models')
    clear models
    
end

% AuReMe
modelDir = 'data/models/aureme';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_no_medium_no_biomass.mat'));
    load(workspace)
    for j=1:numel(models)
        models{j} = addMetFormulae(models{j}, formulaTab);
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    save(workspace, 'models')
    clear models
end