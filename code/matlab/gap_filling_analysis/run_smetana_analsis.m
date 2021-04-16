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
system(['smetana -d ' fullfile(topDir, smetanaDir, 'models', '*.xml') ' -o ' ...
    fullfile(topDir, smetanaDir, 'results', [habitat '_' experiment '_'])])

%% evaulate output
smetanaRes = readtable(fullfile(topDir, smetanaDir, 'results',...
    [habitat '_' experiment '_detailed.tsv']),...
    'ReadVariableNames', true, 'ReadRowNames', false,...
    'Delimiter', '\t', 'FileType', 'text');

% translate metabolite IDs to names
smetanaRes.compound = strtok(smetanaRes.compound, '_');
smetanaRes.cpdName = translateIDs(smetanaRes.compound,'met',[],'MNXref', 'NAME');

otus = unique([smetanaRes.donor; smetanaRes.receiver]);
import_ID = ...
    cellfun(@(x)strjoin(smetanaRes.compound(ismember(smetanaRes.receiver,x)),...
    ','), otus, 'un', 0);
import_NAME = ...
    cellfun(@(x)strjoin(smetanaRes.cpdName(ismember(smetanaRes.receiver,x)),...
    ','), otus, 'un', 0);
export_ID = ...
    cellfun(@(x)strjoin(smetanaRes.compound(ismember(smetanaRes.donor,x)),...
    ','), otus, 'un', 0);
export_NAME = ...
    cellfun(@(x)strjoin(smetanaRes.cpdName(ismember(smetanaRes.donor,x)),...
    ','), otus, 'un', 0);

% create name to BRITE translation
[exchange_name, ia] = unique(smetanaRes.cpdName);
exchange_brite = map2KEGGBrite(smetanaRes.compound(ia), briteFile);
exchange_brite(cellfun(@isempty,exchange_brite)) = {{'Other'}};
exchange_brite = cellfun(@(x)x(1),exchange_brite,'un',0);
exchange_brite = vertcat(exchange_brite{:});

writetable(cell2table([import_ID, import_NAME, export_ID, export_NAME],...
    'VariableNames', {'import_ID', 'import_NAME', 'export_ID', 'export_NAME'},...
    'RowNames', otus), fullfile(topDir, smetanaDir, 'results',...
    ['exchange_' habitat '_' experiment '.txt']), 'WriteVariableNames', true,...
    'WriteRowNames', true, 'Delimiter', '\t');

writetable(cell2table([exchange_brite, exchange_name],...
    'VariableNames', {'BRITE', 'NAME'}), fullfile(topDir, smetanaDir, 'results',...
    ['brite_dict_' habitat '_' experiment '.txt']), 'WriteVariableNames', true,...
    'Delimiter', '\t');


