% Run all evaluation functions on the CarveME reconstructions for all three
% habitats

fprintf('################# Starting Evaluation of CarveMe reconstructions #################\n')
habitats = {'Root', 'Leaf', 'Soil'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n\n', habitats{i})

    
    %% Load the models from workspace
    modelWorkspace = fullfile('data/models/carveme',...
        habitats{i}, strcat(habitats{i}, '_models.mat'));
    if exist(modelWorkspace,'file')
        fprintf('\nWorkspace exists for %s, loading...\n', habitats{i})
        load(modelWorkspace)
    else
        error('No workspace found for %s', habitats{i})
    end
    
    %% Evaluation:
    plotDir = fullfile('figures/CarveMe_draft');
    recMethod = strcat(habitats{i}, '-CarveMe');
    evaluationWorkspace = fullfile('data/models/carveme',...
        habitats{i}, strcat(habitats{i},'_evaluation.mat'));
    evaluateModels(models, recMethod, plotDir, evaluationWorkspace, blackList)
    
end









