% Try fastGapFill

%% Analyze using gapFind
modelPath = '/stud/wendering';
filename = 'Soil522.mat';
load(fullfile(modelPath, filename))

%{
outputMets = detectDeadEnds(model);

DeadEnds = model.mets(outputMets);

[rxnList, rxnFormulaList] = findRxnsFromMets(model, DeadEnds);

model.lb(ismember(model.rxns, rxnList))
model.ub(ismember(model.rxns, rxnList))

[allGaps, rootGaps, downstreamGaps] = gapFind(model, 'true');

BlockedReactions = findBlockedReaction(model);
%}


%% Prepare Stats

Stats = {
    'Model name',...
    'Size S (original model)',...
    'Number of compartments',...
    'List of compartments',...
    'Number of blocked reactions',...
    'Number of solvable blocked reactions',...
    'Size S (flux consistent)',...
    'Size SUX (including solvable blocked reactions)',...
    'Number of added reactions (all)',...
    'Number of added metabolic reactions',...
    'Number of added transport reactions',...
    'Number of added exchange reactions',...
    'Time preprocessing',...
    'Time fastGapFill'};


% Make exchange reactions reversible
EX = strncmp('EX', model.rxns, length('EX'));
model.lb(EX) = -1000;
model.ub(EX) = 1000;
clear EX

% Fill the table with entries
cnt = 1;

% filename
Stats{cnt, 2} = filename; cnt = cnt + 1;

% size of stoichiometric matrix
[a,b] = size(model.S);
Stats{cnt, 2} = strcat(num2str(a),'x',num2str(b)); cnt = cnt + 1;

% number of compartments
[~,rem] = strtok(model.mets, '\[');
rem = unique(rem);
Stats{cnt, 2} = num2str(length(rem)); cnt = cnt + 1;

% list of compartment ids
rem = strjoin(rem, ',');
Stats{cnt, 2} = rem; cnt = cnt + 1;
clear rem

%% Prepare fastGapFill
% 1. flux-consistent model
% 2. SUX matrix: S + universal database (U) placed in all compartments with transport
% reactions for each metabolite (X)

compartments = model.comps;
epsilon = 1E-4;
database = '/net/calc1/srv/wendering/kegg/ligand/reaction/reaction.lst';
dictionary = translateIDs(strtok(model.mets,'\['), 'met', 'ModelSEED', 'KEGG');
dictionaryFile = '/stud/wendering/Masterthesis/DATA/tables/ModelSEED_KEGG_Translation.csv';
writetable(table(strtok(model.mets,'\['), dictionary), dictionaryFile,...
    'Delimiter', '\t', 'WriteVariableNames', false)
blackList = {''};

tic;
[consistModel, consistMatricesSUX, BlockedRxns] = prepareFastGapFill(model,...
    compartments, epsilon, database, dictionaryFile, blackList);
tpre = toc;

% number of all reactions
Stats{cnt, 2} = num2str(length(BlockedRxn.allRxns)); cnt = cnt + 1;

% number of unblocked reactions
Stats{cnt, 2} = num2str(length(BlockedRxns.solvableRxns)); cnt = cnt + 1;

% size of augmented stoichiometrix matrix
[a,b] = size(consistModel.S);
Stats{cnt, 2} = strcat(num2str(a),'x',num2str(b)); cnt = cnt + 1;

% size of SUX matrix
[a,b] = size(consistMatricesSUX.S);
Stats{cnt, 2} = strcat(num2str(a),'x',num2str(b));cnt = cnt + 1;
clear a b
%% Perform fastGapFill
epsilon = 1E-4;

weights.MetabolicRxns = 0.1;
weights.ExchangeRxns = 0.5;
weights.TransportRxns = 10;

tic;
[AddedRxns] = fastGapFill(consistMatricesSUX, epsilon, weights);
tgap = toc;

% number of added reactions
Stats{cnt, 2} = num2str(length(AddedRxns.rxns)); cnt = cnt + 1;

%% Postprocessing

IdentifyPW = 0;
[AddedRxnsExtended] = postProcessGapFillSolutions(AddedRxns, model, BlockedRxns);
clear AddedRxns

% number of added metabolic reactions
Stats{cnt, 2} = num2str(AddedRxnsExtended.Stats.metabolicSol); cnt = cnt + 1;

% number of added transport reaction
Stats{cnt, 2} = num2str(AddedRxnsExtended.Stats.transportSol); cnt = cnt + 1;

% number of added exchange reactions
Stats{cnt, 2} = num2str(AddedRxnsExtended.Stats.exchangeSol); cnt = cnt + 1;

% time required for preparation
Stats{cnt, 2} = num2str(tpre); cnt = cnt + 1;

% time required for gap filling
Stats{cnt, 2} = num2str(tgap); cnt = cnt + 1;

% reaction list
rxnList = {};
col = 1;
rxnList{1, col} = filename; 
rxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxns; col = col + 1;

rxnList{1, col} = filename;
rxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.rxnFormula; col = col + 1;

rxnList{1,col}=filename;
rxnList(2:length(AddedRxnsExtended.rxns)+1,col) = AddedRxnsExtended.subSystem; col = col + 1;

clear AddedRxnsExtended tgap tpre BlockedRxns cnt col consistM*


















    
    
    