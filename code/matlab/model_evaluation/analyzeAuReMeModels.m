% Run all evaluation functions on the AuReMe reconstructions for all three
% habitats
fprintf('################# Starting Evaluation of AuReMe reconstructions #################\n')
habitats = {'Soil', 'Root', 'Leaf'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n', habitats{i})
    %% Load the workspace that contains all models for a given habitat
    % and threshold
    AuReMeModelDir = fullfile('data/models/aureme', habitats{i});
    modelWorkspace = fullfile(AuReMeModelDir, strcat(habitats{i}, '_models.mat'));
    
    if exist(modelWorkspace,'file')
        fprintf('\nWorkspace exists for %s, loading...\n', habitats{i})
        load(modelWorkspace)
    else
        error('No workspace found for %s', habitats{i})
    end
    
    %% Evaluation
    plotDir = 'figures/AuReMe_draft';
    recMethod = strcat(habitats{i}, '-AuReMe');
    evaluationWorkspace = fullfile(AuReMeModelDir,...
        strcat(habitats{i}, '_evaluation.mat'));
    evaluateModels(models, recMethod, plotDir, evaluationWorkspace, blackList)
    
end

