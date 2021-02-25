% Model-generation using the RAVEN toolbox

% Check if installation is correct:

% cd('/stud/wendering/bin/CentOS/RAVEN/software/')
% checkInstallation

% Parameters for model reconstruction

%genomePath  = '/net/calc1/srv/wendering/genomes/ORFs-aminoacid';
path(path, pathdef)
genomePath  = '/srv/wendering/genomes/ORFs-aminoacid';
habitats = {'Leaf', 'Root','Soil'};
dataDir = '/stud/wendering/Masterthesis/DATA/RAVEN_dataDir/selfmade/prok90_kegg90';
outDir = '/srv/wendering/RAVEN-models/outDir';
taxPath = '/stud/wendering/Masterthesis/DATA/genomes';
getModelFromID = 0;
models = cell(numel(habitats),1);
gapFillDir = '/stud/wendering/Masterthesis/DATA/RAVEN_dataDir/template-models';
%modelDir = '/stud/wendering/Masterthesis/DATA/models_RAVEN';
modelDir = '/srv/wendering/RAVEN-models';
%writeToFile = 1;
writeToFile = 2;
p = parpool(8);

tic
for i=1:numel(habitats)
    models{i} = getMultipleModels_par(genomePath, habitats{i}, dataDir, outDir, taxPath, [], getModelFromID, writeToFile, modelDir);
end
save(fullfile(modelDir, 'models_strict_similarity.mat'), models);
toc




%{

%%% Add fillGaps in getMulitpleModels (and addSpontaneous?)

try
    templateModelNames = dir(fullfile(gapFillDir,'*.xml'));
    templateModelNames = {templateModelNames.name};
    templateModels = cell(numel(templateModelNames),1);
catch
    error('The template models for gap filling could not be localized.')
end
for m=1:numel(templateModelNames)
    fprintf('Reading template model %s\n', templateModelNames{m})
    templateModels{m} = importModel(fullfile(gapFillDir,...
        templateModelNames{m}));
end
    


    
    
    
    


% Create a model for Leaf1 using the MetaCyc database

fastaFile = '/stud/wendering/Masterthesis/DATA/genomes/Leaf1.faa';

% BlastP -> WORKS
Leaf1_MetaCyc_BlastP = getMetaCycModelForOrganism('LEAF1',...
    fastaFile,...
    false, false, false, 100, 45, false);


% diamond -> WORKS
Leaf1_MetaCyc_diamond = getMetaCycModelForOrganism('LEAF1',...
    fastaFile);

% Create a model for Leaf1 using the KEGG database -> WORKS
RAVEN_kegg_model_dir = '/stud/wendering/bin/CentOS/RAVEN/external/kegg'
kegg_model_v87 = fullfile(RAVEN_kegg_model_dir,'kegg_model_RAVEN')
kegg_model_v91 = fullfile(RAVEN_kegg_model_dir, 'self-trained-kegg-model')
% Copy self-trained kegg model into the search path of the reconstruction
% function
unix(['cp ', fullfile(kegg_model_v91, '*'), ' ', RAVEN_kegg_model_dir]);

dataDir = '~/Masterthesis/DATA/RAVEN_dataDir/selfmade/prok90_kegg90';
outDir = '~/Masterthesis/DATA/RAVEN_outDir/Leaf1_';

Leaf1_KEGG_self_trained = getKEGGModelForOrganism('Leaf1',...
    fastaFile, dataDir, [], false, false, false, false,...
    10^-30,0.8,0.3,-1);

outDir = '~/Masterthesis/DATA/RAVEN_outDir/Leaf1_noFasta';

Leaf1_KEGG_self_trained_noFasta = getKEGGModelForOrganism('mts',...
    [], dataDir, outDir, false, false, false, false,...
    10^-30,0.8,0.3,-1);




% pre-trained HMM with procaroytic species and 100 % identity
% Copy pre-trained kegg model into the search path of the reconstruction
% function
unix(['cp ', fullfile(kegg_model_v87, '/*', ' ', RAVEN_kegg_model_dir]));

dataDir = '~/Masterthesis/DATA/RAVEN_dataDir/pre-trained/prok100_kegg87';

% pre-trained prokaryotic model -> WORKS
Leaf1_KEGG_pre_trained = getKEGGModelForOrganism('LEAF1',...
    fastaFile, dataDir, outDir, false, false, false, false,...
    10^-50, 0.8, 0.3, -1);










%}





