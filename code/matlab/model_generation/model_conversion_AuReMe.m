% Convert all models to a uniform standard
options; clear
%% AuReMe
% Tables
tablesDir = 'data/tables';
ecTransFile = fullfile(tablesDir, 'corrected-EC-numbers.csv');
taxonomyFile = fullfile(tablesDir, 'taxonomy-units.csv');
formulaeFile = fullfile(tablesDir, 'MNXref', 'MNXref_MET_FORMULAE.csv');

formulaTab = readtable(formulaeFile, 'ReadVariableNames', true, 'Delimiter', '\t');
taxonomyTab = readtable(taxonomyFile, 'ReadVariableNames', false);
ecTranslationTable = readtable(ecTransFile, 'ReadVariableNames', false);
translationDB = loadTranslationDB;


fprintf('################# Converting AuReMe draft models #################\n')
habitats = {'Leaf', 'Root', 'Soil'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n', habitats{i})
    
    % Model directory
    AuReMeModelDir = fullfile('data/models/aureme',habitats{i});
    ModelWorkspace = fullfile('data/models/aureme',...
    habitats{i}, strcat(habitats{i}, '_models.mat'));
    DraftWorkspace = fullfile('data/models/aureme',...
        habitats{i},strcat(habitats{i},'_models_draft.mat'));
    
    AuReMeModels = dir(fullfile(AuReMeModelDir, '*.sbml'));
    AuReMeModels = {AuReMeModels.name};
    names = cellfun(@(x)regexp(x, strcat(habitats{i},'[0-9]*[^\.]*'), 'match'), AuReMeModels);
    
    n = numel(AuReMeModels);
    
    % Reading from file of draft models are not saved in a workspace yet
    if ~isfile(DraftWorkspace)
        for j=1:n
            fprintf('Reading file #%d/%d: %s\n', j, n, AuReMeModels{j})
            eval(strcat(names{j},'=','readSBML(fullfile(AuReMeModelDir,AuReMeModels{j}), 1000);'))
        end
        
        
        % Collect the single models into a cell array 'models' and delete the
        % individual models
        eval(strcat('models = {', strjoin(names, ';'), '};'))
        eval(string(strcat('clear', {' '}, habitats{i}, '*')))
        
        % make the stoichiometric matrices sparse to be able to save the models
        % properly
        for j=1:n
            models{j}.S = sparse(models{j}.S);
        end
        
        save(DraftWorkspace, 'models');
        
    elseif ~isfile(ModelWorkspace)
        load(DraftWorkspace)
        n = numel(models);
    else 
        load(ModelWorkspace)
        n = numel(models);
    end
    
    if ~isfield(models{end}, 'EC')
        % re-order taxonomy names according to the order of models
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
        
        % Add E.C. numbers
        fprintf('\nAdding EC numbers\n')
        rxnTab = translationDB.rxnTab;
        parfor j=1:n
            % Retrieve ec numbers from ModelSEED reaction database
            ec = translateIDs(models{j}.rxns, 'rxn', rxnTab, 'MetaCyc', 'EC');
            % Correct the EC numbers as they could be outdated
            models{j}.EC =  correctEC(ec, ecTranslationTable);
            % add a description
            models{j}.description = strcat(strtok(tax_names{j}, '_'), '_AuReMe_model');
            % add a formula
            models{j} = addMetFormualae(models{j}, formulaTab);
        end
        clear metTab ec
        save(ModelWorkspace, 'models');
        fprintf('Done.\n')
    end
    if ~any(contains(models{end}.rxns, 'MNXR'))
	% convert the model to the the standard format
	fprintf('\nConverting to standard format and translating to MNXref namespace\n')
    	for j=1:n
            models{j} = convertAuReMeModel(models{j}, true, translationDB, false);
    	end
    else
	fprintf('Already converted')
    end
    save(ModelWorkspace, 'models');
end

