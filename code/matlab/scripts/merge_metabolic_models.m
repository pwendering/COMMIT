% merge draft metabolic models from different approaches
c = parcluster;
c.NumWorkers = 4;
delete(gcp('nocreate'))
P = parpool(c);
habitats = {'Soil', 'Leaf', 'Root'};
methods = {'CarveMe', 'KBase', 'AuReMe', 'RAVEN'};
dataDir = '/stud/wendering/Masterthesis/DATA';
disp('-------------------------------------------------------------------')
disp('START')
disp('-------------------------------------------------------------------')
disp('loading the universal database')
load('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat')
disp('-------------------------------------------------------------------')
for i=1:numel(habitats)
    disp(habitats{i})
    disp('------------------------------')
    disp('Loading and collecting models from the different approaches...')
    for j=1:numel(methods)
        if isequal(methods{j}, 'RAVEN')
            workspace = fullfile(dataDir, strcat('models_', methods{j}),...
                'HMMer10E-50', strcat(habitats{i}, '_models_COBRA_GPR'));
        else
            workspace = fullfile(dataDir, strcat('models_', methods{j}),...
                habitats{i}, strcat(habitats{i}, '_models_metFormulas'));
        end
        load(workspace);
        eval([strcat(habitats{i}, '_', methods{j}), '= models;']);
    end
    disp('------------------------------')
    merged_models = {};
    for j=1:numel(models)
        id = strtok(models{j}.id, '_');
        fprintf('Model #%d (%s)\n', j, id)
        fprintf('Collecting %d models...\n', numel(methods))
        % create models variable for every OTU in each habitat
        models_to_merge = {};
        for k=1:numel(methods)
            eval(['models_to_merge = vertcat(models_to_merge,',...
                strcat(habitats{i}, '_', methods{k}), '{j});'])
        end
        disp('------------------------------')
        merged_models{j} = mergeModels(models_to_merge, dbModel_MNXref_balanced);
        wo_del = sum(cellfun(@(x)numel(x.rxns), models_to_merge));
        w_del = numel(merged_models{j}.rxns);
        fprintf('Number of deleted (merged) reactions:\t%d\n', wo_del-w_del)
        disp('------------------------------')
    end
    
    workspace = fullfile(dataDir, 'Consensus_models',...
        strcat(habitats{i}, '_consensus_models'));
    models = models_to_merge;
    disp('saving workspace')
    save(workspace, 'merged_models')
    clear models
    disp('-------------------------------------------------------------------') 
end
        
     
