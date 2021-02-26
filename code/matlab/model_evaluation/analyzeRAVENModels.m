% Run all evaluation functions on the RAVEN reconstructions for all three
% habitats

fprintf('################# Starting Evaluation of RAVEN reconstructions #################\n')
habitats = {'Leaf', 'Root', 'Soil'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n', habitats{i})
    
    %% Load the workspace that contains all models for a given habitat
    RAVENModelDir = fullfile(topDir, 'data/models/raven/HMMer10E-50');
    modelWorkspace = fullfile(RAVENModelDir,...
        strcat(habitats{i}, '_converted.mat'));
    
    if exist(modelWorkspace,'file')
        fprintf('\nWorkspace exists for %s, loading...\n', habitats{i})
        load(modelWorkspace)
    else
        error('No workspace found for %s', habitats{i})
    end
    
    %% Evaluation
    plotDir = fullfile(topDir, 'figures/RAVEN_draft');
    recMethod = strcat(habitats{i},'-HMMer10E-50-RAVEN');
    evaluationWorkspace = fullfile(RAVENModelDir,...
        strcat(habitats{i}, '_evaluation.mat'));
    evaluateModels(models, recMethod, plotDir, evaluationWorkspace, blackList)
    
    
    
end
