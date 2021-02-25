% Convert all models to a uniform standard

%% RAVEN
% Tables
tablesDir = '/stud/wendering/Masterthesis/DATA/tables';
ecTransFile = fullfile(tablesDir, 'corrected-EC-numbers.csv');
coFactorFile = fullfile(tablesDir, 'cofactors_from_KEGG.csv');
taxonomyFile = fullfile(tablesDir, 'taxonomy-units.csv');
metTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-met-translation-table.csv');
rxnTransFile = fullfile(tablesDir, 'MNXref', 'MNXref-rxn-translation-table.csv');
RAVENModelDir = '/srv/wendering/RAVEN-models';
formulaeFile = fullfile(tablesDir, 'MNXref', 'MNXref_MET_FORMULAE.csv');

formulaTab = readtable(formulaeFile, 'ReadVariableNames', true, 'Delimiter', '\t');
ecTranslationTable = readtable(ecTransFile, 'ReadVariableNames', false);
coFactorsTab = readtable(coFactorFile, 'Delimiter', '\t', 'ReadVariableNames', false);
taxonomyTab = readtable(taxonomyFile, 'ReadVariableNames', false);
metTransTab = readtable(metTransFile, 'ReadVariableNames', true);
rxnTransTab = readtable(rxnTransFile, 'ReadVariableNames', true);
translationDB.metTab = metTransTab;
translationDB.rxnTab = rxnTransTab;

fprintf('################# Converting RAVEN draft models #################\n')
habitats = {'Root', 'Soil'};%, 'Leaf'};
thresholds = {'30'};%{'30', '50'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n', habitats{i})
    for t=1:numel(thresholds)
        fprintf('\n################# HMMER threshold: 10E-%s\n', thresholds{t})
        
        % Load the workspace that contains all models for a given habitat
        % and threshold
        ws = fullfile(RAVENModelDir, strcat('HMMER-threshold_10E-',...
            thresholds{t}), strcat(habitats{i}, '.mat'));
        fprintf('\nLoading workspace from %s\n', ws)
        load(ws);
        
        % Get all objects in the workspace that contain that habitat name
        names =  who('-regexp',habitats{i});
        
        % Store the models in a cell array for better processing and delete
        % the single models
        eval(strcat('models = {', strjoin(names, ';'),'};'))
        eval(string(strcat('clear',{' '},habitats{i},'*')))
        
        % get taxonomy names
        [~,idx] = ismember(names, taxonomyTab.Var1);
        tax_names = taxonomyTab.Var2(idx);
        
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
        
        % convert the model to the the standard format
        fprintf('\nConverting to standard format and translating to MNXref namespace\n')
        parfor j=1:numel(models)
            models{j} = convertRavenModel(models{j}, true, translationDB);
            models{j}.EC = correctEC(models{j}.EC, ecTranslationTable);
            models{j}.description = strcat(strtok(tax_names{j}, '_'), '_RAVEN_model');
            % add a formula
            models{j} = addMetFormualae(models{j}, formulaTab);
        end
        
        save(fullfile(RAVENModelDir, strcat('HMMER-threshold_10E-',...
            thresholds{t}), strcat(habitats{i}, '_converted.mat')), 'models');
    end
end


