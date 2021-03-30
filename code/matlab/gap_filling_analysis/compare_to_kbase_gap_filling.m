% compare gap-filling solutions of ComGapFill with the results from the
% KBase gap-filling App (using auxtropic media)
experiment = 'Bulgarelli';
tableOutFile = ['data/gap-filling/iterative/Soil/KBase/', experiment,...
    '_KBase_comp.csv'];

% list model files 
modelDir = 'data/models/kbase/Soil';
dirContent = dir(fullfile(modelDir,'Soil-models-gf-auxo','*.sbml'));
modelFiles = {dirContent.name};
modelIDs = strtok(modelFiles,'.');

% find added reactions by comparing the reaction sets of the draft model
% with the one of the gap-filled model
disp('Reading KBase models')
gfRxns = cell(1,numel(modelFiles));
for i=1:numel(modelFiles)
    fprintf('Model #%d...', i)
    
    model1 = readCbModel(fullfile(modelDir, 'Soil-models-gf-auxo', modelFiles{i}));
    model2 = readCbModel(fullfile(modelDir, 'Soil-models', modelFiles{i}));
    
    gfRxns{i} = setdiff(model1.rxns,model2.rxns);
    fprintf('done\n')
end

% translate IDs to MNXref namespace
disp('Translating IDs')
gfRxnsMNX = cell(1,numel(modelFiles));
excMets = cell(1,numel(modelFiles));

for i=1:numel(gfRxns)
    rxnIDs = regexp([gfRxns{i}{:}],'rxn\d{5}','match')';
    metIDs = regexp([gfRxns{i}{:}],'cpd\d{5}','match')';
    gfRxnsMNX{i} = translateIDs(rxnIDs, 'rxn',[],'ModelSEED','MNXref',false);
    excMets{i} = translateIDs(metIDs, 'met', [], 'ModelSEED','MNXref', false);
end

clearvars -except gfRxnsMNX gfRxns excMets modelIDs experiment

% find added reactions of gap-filled KBase models using ComGapFill
disp('Find added reactions in ComGapFilled KBase models')
modelFile = ['data/gap-filling/iterative/Soil/KBase/' experiment '.mat'];
load(modelFile,'GF')

gfRxnsComGf = cell(1, numel(GF));
for i=1:numel(GF)
    gfRxnsComGf{i} = GF{i}.rxns(ismember(cellfun(@char,GF{i}.rxnNotes,'un',0),'gf'));
end

% compare reaction sets and numbers of added reactions
disp('Comparing added reactions')
nKBase = zeros(numel(GF),1);
nComGF = zeros(numel(GF),1);
nExcKBase = zeros(numel(GF),1);
nExcComGF = zeros(numel(GF),1);
gfIntersect = cell(numel(GF),1);
excMetIntersect = cell(numel(GF),1);
IDs = repmat({''},1,numel(GF));

for i=1:numel(GF)
    IDs{i} = strtok(GF{i}.id, '_');
    idxMatch = ismember(modelIDs, strtok(GF{i}.id, '_'));
    
    gfIntersect{i} = intersect(strtok(gfRxnsComGf{i},'_'), gfRxnsMNX{idxMatch});
    excMetIntersect{i} = intersect(...
        regexp([gfRxnsComGf{i}{:}],'MNXM\d+','match'), excMets{idxMatch});
    nKBase(i) = numel(gfRxnsMNX{idxMatch});
    nComGF(i) = numel(gfRxnsComGf{i});
    nExcKBase(i) = sum(contains(gfRxns{idxMatch},'EX_'));
    nExcComGF(i) = sum(contains(gfRxnsComGf{i}, 'EX_'));
end


% write results to file
resultTable = array2table([nKBase, nComGF, cellfun(@numel, gfIntersect),...
    nExcKBase, nExcComGF, cellfun(@numel, excMetIntersect)],...
    'VariableNames', {'KBase', 'ComGapFill', 'Intersection',...
    'EX_KBase', 'EX_ComGapFill', 'EX_Intersection'},...
    'RowNames', IDs);
writetable(resultTable, tableOutFile, 'WriteVariableNames', true,...
    'WriteRowNames', true);
