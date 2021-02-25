% add formulae to models

formulaTab = readtable('/stud/wendering/Masterthesis/DATA/tables/MNXref/MNXref_MET_FORMULAE.csv', 'ReadVariableNames', true);
habitats = {'Soil', 'Leaf', 'Root'};

% KBase models
modelDir = '/stud/wendering/Masterthesis/DATA/models_KBase';
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
thresholds = {'50'};%, '30'};

for t = thresholds
    modelDir = fullfile('/stud/wendering/Masterthesis/DATA/models_RAVEN', ...
        strcat('HMMer10E-', t{:}));
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
    
end

% CarveMe
modelDir = '/stud/wendering/Masterthesis/DATA/models_CarveMe';
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
modelDir = '/stud/wendering/Masterthesis/DATA/models_AuReMe';
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