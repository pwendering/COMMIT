% Options for gap filling

%% Gap filling

%~~~~~~~~~~~~~ General settings ~~~~~~~~~~~~~%

% handling of warnings
warning('off', 'all')

% number of parallel workers
ncpu=8;

% habitat
habitat = 'Leaf';

% experiments to be analyzed
experiments = {'leaf_natural_sites', 'IPL'};

% medium file
mediumFile = '/stud/wendering/Masterthesis/DATA/media/minimal-medium.mat';

% directory for auxotrophy media (files must be of form 'model id'.tsv)
mediaDir = fullfile('/stud/wendering/Masterthesis/DATA/media', habitat, 'auxo-media');

% top directory for OTU abundance (subdirectories must be habitat >
% experiment > otutab.txt
otuDir = '/stud/wendering/Masterthesis/DATA/Abundances/';

% threshold for the biomass reaction
epsilon = 1E-4;

% number of iterations for the iterative approach
iterations = 100;

% print level
verbose = true;

% provide an order for iterative gap filling
order = [];

% whether or not sink reactions should be included in the gap-filling
% objective
include_sink = false;

% add path for IBM CPLEX solver
addpath('/stud/wendering/bin/Fedora/CPLEX-Optimizer/cplex/matlab/x86-64_linux')

%~~~~~~~~~~~~~ Model workspace and output directory ~~~~~~~~~~~~~%

% location where draft models are stored
modelDir = '/stud/wendering/Masterthesis/DATA/Consensus_models';
% modelFile = fullfile(modelDir, [habitat, '_consensus_models_noCarveMe_biomass']);
modelFile = fullfile(modelDir, [habitat, '_consensus_models_biomass']);
% modelDir = '/stud/wendering/Masterthesis/DATA/models_KBase';
% modelFile = fullfile(modelDir, habitat, [habitat, '_models_biomass']);

tmp_spec = '_all';

% output directory for the gap-filled models
% outDir = fullfile('/stud/wendering/Masterthesis/DATA/Gap-filling/iterative', habitat, '_KBase');
% outDir = fullfile('/stud/wendering/Masterthesis/DATA/Gap-filling/iterative', habitat, '_no_CarveMe');
outDir = fullfile('/stud/wendering/Masterthesis/DATA/Gap-filling/iterative', habitat, 'all');


%~~~~~~~~~~~~~ Gap-filling resources ~~~~~~~~~~~~~%

% file for the database
dbFile = '/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat';

% path to sequence similarity workspace
seq_sim_workspace = '/stud/wendering/Masterthesis/DATA/Gap-filling/sequence-similarity/sequence_similarity.mat';


%~~~~~~~~~~~~~ Reaction-type-specific weights ~~~~~~~~~~~~~%

% reactions with sequence evidence (E >= 10E-6) 
weights.sequence = 25;
% transport reactions
weights.transport = 100;
% transport reactions operating on highly-permeable metabolites
weights.permeable = 50;
% general weight for metabolic reactions
weights.metabolic = 50;
% weight for making a reaction reversible
weights.reverse = 25;
% weight for exchange reactions from the database
weights.exchange = 100000;
% weight for allowed exchange / uptake reactions (from previous solutions)
weights.uptake = 1;
% weight for sink reactions if they should be determined in the gap
% filling program
weights.sink = 0;
% weights for reactions already contained in the model
weights.model = 0;

