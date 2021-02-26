% merge draft metabolic models from different approaches
options
clearvars -except topDir dbFile ncpu

% set up parallel pool
c = parcluster;
c.NumWorkers = ncpu;
delete(gcp('nocreate'))
P = parpool(c);

habitats = {'Soil', 'Leaf', 'Root'};
methods = {'kbase', 'aureme', 'raven'};
modelTopDir = fullfile(topDir, 'data', 'models');

disp('-------------------------------------------------------------------')
disp('START')
disp('-------------------------------------------------------------------')

disp('loading the universal database')
load(dbFile)

disp('-------------------------------------------------------------------')
for i=1:numel(habitats)
    disp(habitats{i})
    disp('------------------------------')
    disp('Loading and collecting models from the different approaches...')
    for j=1:numel(methods)
        if isequal(methods{j}, 'raven')
            workspace = fullfile(modelTopDir, methods{j},...
                'HMMer10E-50', strcat(habitats{i}, '_models_COBRA_GPR'));
        else
            workspace = fullfile(modelTopDir, methods{j},...
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
        
        % run merging function
        merged_models{j} = mergeModels(models_to_merge, dbModel_MNXref_balanced);
        wo_del = sum(cellfun(@(x)numel(x.rxns), models_to_merge));
        w_del = numel(merged_models{j}.rxns);
        fprintf('Number of deleted (merged) reactions:\t%d\n', wo_del-w_del)
        disp('------------------------------')
        
    end
    
    workspace = fullfile(modelTopDir, 'consensus',...
        strcat(habitats{i}, '_consensus_models_noCarveMe'));
    models = models_to_merge;
    disp('saving workspace')
    save(workspace, 'merged_models')
    clear models
    disp('-------------------------------------------------------------------') 
end
clear topDir
     
