% Convert all models to a uniform standard

%% KBase
% Tables
tablesDir = '/stud/wendering/Masterthesis/DATA/tables';
ModelSEEDRxnEC = fullfile(tablesDir, 'ModelSEED_RXN_EC.csv');
ecTransFile = fullfile(tablesDir, 'corrected-EC-numbers.csv');
metTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-met-translation-table.csv');
rxnTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-rxn-translation-table.csv');
formulaeFile = fullfile(tablesDir, 'MNXref', 'MNXref_MET_FORMULAE.csv');

formulaTab = readtable(formulaeFile, 'ReadVariableNames', true, 'Delimiter', '\t');
ecRxnTable = readtable(ModelSEEDRxnEC, 'ReadVariableNames', false);
ecTranslationTable = readtable(ecTransFile, 'ReadVariableNames', false);
metTransTab = readtable(metTransFile, 'ReadVariableNames', true);
rxnTransTab = readtable(rxnTransFile, 'ReadVariableNames', true);
translationDB.metTab = metTransTab;
translationDB.rxnTab = rxnTransTab;


fprintf('################# Converting KBase draft models #################\n')
habitats = {'Leaf', 'Root', 'Soil'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n', habitats{i})
    
    % Model directory
    KBaseModelDir = fullfile('/srv/wendering/KBase-models',habitats{i});
    ModelWorkspace = fullfile('/stud/wendering/Masterthesis/DATA/models_KBase',...
    habitats{i},strcat(habitats{i},'_models.mat'));
    KBaseModels = dir(fullfile(KBaseModelDir, '*.sbml'));
    KBaseModels = {KBaseModels.name};
    names = cellfun(@(x)regexp(x, strcat(habitats{i},'[0-9]*[^\.]*'), 'match'),KBaseModels);
    
    n = numel(KBaseModels);
    
    % Reading from file
    for j=1:n
        fprintf('Reading file #%d/%d: %s\n', j, n, KBaseModels{j})
        eval(strcat(names{j},'=','readSBML(fullfile(KBaseModelDir,KBaseModels{j}), 1000);'))
    end
    
    % Collect the single models into a cell array 'models' and delete the
    % individual models
    eval(strcat('models = {', strjoin(names, ';'),'};'))
    eval(string(strcat('clear',{' '},habitats{i},'*')))
    
    % make the stoichiometric matrices sparse to be able to save the models
    % properly
    for j=1:n
        models{j}.S = sparse(models{j}.S);
    end
    
    save(ModelWorkspace, 'models');
    
    % Add E.C. numbers
    fprintf('\nAdding EC numbers\n')
    parfor j=1:n
        %fprintf('\mModel %d/%d:\n', num2str(j), num2str(n))
        % Retrieve ec numbers from ModelSEED reaction database
        ec = ecNumbersFromModelSEED(models{j}.rxns, ecRxnTable);
        % Correct the EC numbers as they could be outdated
        models{j}.EC =  correctEC(ec, ecTranslationTable);
    end
    save(ModelWorkspace, 'models');
    
    clear ec
    fprintf('Done.\n')
    
    % convert the model to the the standard format
    fprintf('\nConverting to standard format and translating to MNXref namespace\n')
    parfor j=1:n
        %fprintf('\mModel %d/%d:\n', num2str(j), num2str(n))
        models{j} = convertKBaseModel(models{j}, true, translationDB);
        % add a formula
        models{j} = addMetFormualae(models{j}, formulaTab);
    end
    
    save(ModelWorkspace, 'models');
end

