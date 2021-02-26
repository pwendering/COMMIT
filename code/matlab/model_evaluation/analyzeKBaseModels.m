% Run all evaluation functions on the KBase reconstructions for all three
% habitats

fprintf('################# Starting Evaluation of KBase reconstructions #################\n')
habitats = {'Soil', 'Root', 'Leaf'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n', habitats{i})
     
   %% Load the models from workspace
    modelWorkspace = fullfile(topDir, 'data/models/kbase',...
        habitats{i},strcat(habitats{i},'_models.mat'));
    if exist(modelWorkspace,'file')
        fprintf('\nWorkspace exists for %s, loading...\n', habitats{i})
        load(modelWorkspace)
    else
        error('No workspace found for %s', habitats{i})
    end
    
    %% Evaluation:
    plotDir = fullfile(topDir, 'figures/KBase_draft');
    recMethod = strcat(habitats{i}, '-Kbase');
    evaluationWorkspace = fullfile(topDir, 'data/models/kbase',...
        habitats{i}, strcat(habitats{i},'_evaluation.mat'));
    evaluateModels(models, recMethod, plotDir, evaluationWorkspace, blackList)
    
end






