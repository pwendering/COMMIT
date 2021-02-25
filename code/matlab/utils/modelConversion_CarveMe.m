% Convert all models to a uniform standard

%% CarveMe
% Tables
tablesDir = '/stud/wendering/Masterthesis/DATA/tables';
ecTransFile = fullfile(tablesDir, 'corrected-EC-numbers.csv');
coFactorFile = fullfile(tablesDir, 'cofactors_from_KEGG.csv');
taxonomyFile = fullfile(tablesDir, 'taxonomy-units.csv');
metTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-met-translation-table.csv');
rxnTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-rxn-translation-table.csv');
formulaeFile = fullfile(tablesDir, 'MNXref', 'MNXref_MET_FORMULAE.csv');

formulaTab = readtable(formulaeFile, 'ReadVariableNames', true, 'Delimiter', '\t');
ecTranslationTable = readtable(ecTransFile, 'ReadVariableNames', false);
coFactorsTab = readtable(coFactorFile, 'Delimiter', '\t', 'ReadVariableNames', false);
taxonomyTab = readtable(taxonomyFile, 'ReadVariableNames', false);
metTransTab = readtable(metTransFile, 'ReadVariableNames', true);
rxnTransTab = readtable(rxnTransFile, 'ReadVariableNames', true);
translationDB.metTab = metTransTab;
translationDB.rxnTab = rxnTransTab;

fprintf('################# Converting CarveMe draft models #################\n')
habitats = {'Root', 'Soil', 'Leaf'};
for i=1:numel(habitats)
    
    fprintf('\n################# %s\n\n', habitats{i})
    
    CarveMeModelDir = fullfile('/stud/wendering/Masterthesis/DATA/models_CarveMe',habitats{i});
    DraftWorkspace = fullfile('/stud/wendering/Masterthesis/DATA/models_CarveMe',...
        habitats{i},strcat(habitats{i},'_models_draft.mat'));
    ModelWorkspace = fullfile('/stud/wendering/Masterthesis/DATA/models_CarveMe',...
        habitats{i},strcat(habitats{i},'_models.mat'));
    CarveMeModels = dir(fullfile(CarveMeModelDir, '*.xml'));
    CarveMeModels = {CarveMeModels.name};
    names = cellfun(@(x)regexp(x, strcat(habitats{i},'[0-9]*[^\.]*'), 'match'),CarveMeModels);
    
    n = numel(CarveMeModels);
    
    for j=1:n
        fprintf('Reading file #%d/%d: %s\n', j, n, CarveMeModels{j})
        eval(strcat(names{j},'=','readCbModel(fullfile(CarveMeModelDir,CarveMeModels{j}));'))
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
    
    save(DraftWorkspace, 'models')
    
    % get taxonomy names
    [~,idx] = ismember(names, taxonomyTab.Var1);
    tax_names = taxonomyTab.Var2(idx);
    clear idx
    % Add numbers to duplicate names
    fprintf('\nProcessing taxonomy names\n')
    for j=1:numel(models)
        count = 1;
        name = strcat(tax_names{j},'_', num2str(count));
        while ismember(name, tax_names)
            count = count + 1;
            name = regexpi(name, '[a-z]*', 'match');
            name = strcat(name{:}, '_',  num2str(count));
        end
        tax_names{j} = name;
    end
    clear count name
    fprintf('Done.\n')
    fprintf('\nConverting to standard format and translating to MNXref namespace\n')
    
    parfor j=1:numel(models)
        models{j} = convertCarveMeModel(models{j}, true, translationDB);
        models{j}.EC = correctEC(models{j}.EC, ecTranslationTable);
        models{j}.description = strcat(strtok(tax_names{j}, '_'), '_CarveMe_model');
        % add a formula
        models{j} = addMetFormualae(models{j}, formulaTab);
    end
    
    save(ModelWorkspace, 'models')
    
end