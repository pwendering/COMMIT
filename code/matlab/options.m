% Options for gap filling

%% Gap filling

%~~~~~~~~~~~~~ General settings ~~~~~~~~~~~~~%

% handling of warnings
warning('off', 'all')

% number of parallel workers
ncpu = 8;

% habitat
habitat = 'Leaf';

% experiments to be analyzed
experiments = {'leaf_natural_sites', 'IPL'};

% path to ComGapFil top directory
global topDir
topDir = '/stud/wendering/ComGapFill';

% medium file
mediumFile = fullfile(topDir, 'data/media/minimal-medium.mat');

% directory for auxotrophy media (files must be of form 'model id'.tsv)
mediaDir = fullfile(topDir, 'data/media', habitat, 'auxo-media');

% top directory for OTU abundance (subdirectories must be habitat >
% experiment > otutab.txt
otuDir = fullfile(topDir, 'data/otu_abundances/');

% threshold for the biomass reaction
epsilon = 1E-4;

% number of iterations for the iterative approach
iterations = 100;

% print level
verbose = true;

% provide an order for iterative gap filling (optional)
order = [];

% whether or not sink reactions should be included in the gap-filling
% objective
include_sink = false;

% set COBRA solver for linear optimization problems
changeCobraSolver('ibm_cplex','LP',0);

%~~~~~~~~~~~~~ Model workspace and output directory ~~~~~~~~~~~~~%

% location where draft models are stored
modelDir = fullfile(topDir, 'data/models/consensus/');
% modelFile = fullfile(modelDir, [habitat, '_consensus_models_noCarveMe_biomass']);
modelFile = fullfile(modelDir, [habitat, '_consensus_models_biomass']);
% modelDir = fullfile(topDir, 'data/models/kbase');
% modelFile = fullfile(modelDir, habitat, [habitat, '_models_biomass']);

tmp_spec = '_all';
sub_dir = strtok(tmp_spec,'_');

% output directory for the gap-filled models
% outDir = fullfile(topDir, 'data/gap-filling/iterative', habitat, sub_dir);
% outDir = fullfile(topDir, 'data/gap-filling/iterative', habitat, sub_dir);
outDir = fullfile(topDir, 'data/gap-filling/iterative', habitat, sub_dir);


%~~~~~~~~~~~~~ Gap-filling resources ~~~~~~~~~~~~~%

% file for the database
dbFile = fullfile(topDir, 'data/gap-filling/database/Universal-model-MNXref-balanced.mat');

% path to sequence similarity workspace
seq_sim_workspace = fullfile(topDir, 'data/gap-filling/sequence-similarity/sequence_similarity.mat');


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
% weight for sink reactions if they should be determined in the gap-filling program
weights.sink = 0;
% weights for reactions already contained in the model
weights.model = 0;

%~~~~~~~~~~~~~ metabolite classification resources ~~~~~~~~~~~~~%

% ChEBI ontology DAG workspace
chebiOntologyWS = fullfile(topDir, 'data/metabolite-classification/ontologyGraph.mat');

% tab-separated ontology file
ontologyFile = fullfile(topDir, 'data/metabolite-classification/chebi_ontology.csv');

% brite Hierarchy file
briteFile = fullfile(topDir, 'data/metabolite-classification/briteHierarchy_ext.csv');
