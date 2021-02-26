% Translate gene IDs 
options
clearvars -except topDir
habitats = {'Soil', 'Leaf', 'Root'};

%% AuReMe
modelDir = fullfile(topDir, 'data/models/aureme');
blast_res_dir = fullfile(topDir, 'data/genomes/DFast-annotations-blast-results-evalue-10');

for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        filename = fullfile(blast_res_dir, habitats{i},...
            strcat(model.id, '.mapping'));
        model= translateGeneIDs(model, filename);
        models{j} = model; clear model
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_genes_translated.mat'));
    save(workspace, 'models')
    clear models
end

%% KBase
modelDir = fullfile(topDir, 'data/models/kbase');
blast_res_dir = fullfile(topDir, 'data/genomes/KBase-annotations-blast-results');

for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        id = strtok(model.id, '_');
        filename = fullfile(blast_res_dir, habitats{i},...
            strcat(id, '.mapping'));
        model= translateGeneIDs(model, filename);
        models{j} = model; clear model
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_genes_translated.mat'));
    save(workspace, 'models')
    clear models
end

