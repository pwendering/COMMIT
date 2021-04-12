% Run SMETANA analysis
% Algorithm published by Zelesniak et al. (2015) PNAS
% Python implementation by Daniel Machado
% (https://github.com/cdanielmachado/smetana)
options

%% Input
% model workspace
habitat = 'Soil';
spec = 'all';
experiment = 'Schlaeppi';
modelFile = fullfile('data/gap-filling/iterative/',...
    habitat, spec, experiment);
load(modelFile, 'GF', 'EX');
n = numel(GF);

smetanaDir = 'smetana-analysis';
% check whether model and result directories exist
dirs = {'models', 'results'};
for i=1:numel(dirs)
    if ~exist(fullfile(smetanaDir,dirs{i}),'dir')
        mkdir(fullfile(smetanaDir,dirs{i}))
    end
end

%% write medium to tsv file
EX = setdiff(EX, 'BIOMASS[e]');
EX = strtok(EX,'[');

writetable(cell2table(...
    [repmat({'min','minimal medium'}, numel(EX), 1), EX', EX'],...
    'VariableNames', {'medium', 'description', 'compound', 'name'}),...
    fullfile(smetanaDir, 'media', [strjoin({habitat experiment spec}, '_') '.tsv']),...
    'FileType', 'text', 'Delimiter', '\t')

%% Prepare models
fprintf('\nWriting %d models to .SBML files\n',n)
for i=1:n
    fprintf('Processing model #%d...',i)
    % reduce to minimum
    tmpModel = rmfield(GF{i}, {'subSystems', 'EC',...
        'grRules', 'rxnNotes', 'rxnNames'});
    for j=1:numel(tmpModel.metNames)
        tmpModel.metNames(j) = cellstr(tmpModel.metNames{j});
    end
    
    tmpModel = convertModelToReversible(tmpModel);
    
    
    % remove obsolete exchange reacions
    removeIdx = ~cellfun(@isempty, regexp(tmpModel.rxns, '^EX_.*_r$'));
    
    for j=1:numel(tmpModel.rxns)
        if contains(tmpModel.rxns(j), 'EX_')
            % modify compartment encoding
            tmpModel.rxns(j) = strcat(regexprep(tmpModel.rxns(j),'\[.\]',''));
            tmpModel.rxns{j} = [tmpModel.rxns{j} '_e'];
            % make reversible, unidirectional reaction in opposite
            % direction will be removed
            tmpModel.lb(j) = -1000;           
        end
        
    end
    

    
    % remove export and sink reactions introduced during gap filling
%     removeIdx = removeIdx | contains(tmpModel.rxns, 'export_M') | contains(tmpModel.rxns, 'sink_M');
    
    tmpModel = removeRxns(tmpModel, tmpModel.rxns(removeIdx));
    
    tmpModel.metMetaNetXID = strtok(tmpModel.mets, '[');
    
    % remove blocked reactions
    blocked = findBlockedReaction(tmpModel);
    tmpModel = removeRxns(tmpModel, blocked);
    
    % write model in SBML format
    writeCbModel(tmpModel, 'format', 'sbml',...
        'fileName', fullfile(smetanaDir, 'models', tmpModel.id));
    
end

%% call smetana from command line
% detailed analysis
system(['smetana -d ' fullfile(topDir, smetanaDir, 'models', '*') ' -o ' ...
    fullfile(topDir, smetanaDir, 'results')])

%% evaulate output