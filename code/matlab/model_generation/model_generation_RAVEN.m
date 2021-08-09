% Model-generation using the RAVEN toolbox
options; clear

% Check if installation is correct:

% cd('~/bin/RAVEN/software/')
% checkInstallation

% Parameters for model reconstruction
genomePath  = 'data/genomes';
orfPath = fullfile(genomePath, 'ORFs-aminoacid');
habitats = {'Leaf', 'Root','Soil'};
dataDir = 'data/RAVEN_dataDir/selfmade/prok90_kegg90';
outDir = 'data/RAVEN_outDir';
getModelFromID = 0;
models = cell(numel(habitats),1);
gapFillDir = 'data/RAVEN_dataDir/template-models';
modelDir = 'data/models/raven';

writeToFile = 2;
p = parpool(ncpu);

tic
for i=1:numel(habitats)
    models{i} = getMultipleModels(orfPath, habitats{i}, dataDir, outDir, genomePath, [], getModelFromID, writeToFile, modelDir);
end
save(fullfile(modelDir, 'models_strict_similarity.mat'), models);
toc