% Run all evaluation functions on the RAVEN reconstructions for all three
% habitats

fprintf('################# Starting Evaluation of RAVEN reconstructions #################\n')
habitats = {'Leaf'};%, 'Root', 'Soil'};
thresholds = {'50'};%, '30'};

for i=1:numel(habitats)
    fprintf('\n################# %s\n', habitats{i})
    for t=1:numel(thresholds)
        fprintf('\n################# HMMER threshold: 10E-%s\n', thresholds{t})
        
        %% Load the workspace that contains all models for a given habitat
        % and threshold
        RAVENModelDir = fullfile('/stud/wendering/Masterthesis/DATA/models_RAVEN',...
            strcat('HMMer10E-', thresholds{t}));
        modelWorkspace = fullfile(RAVENModelDir,...
            strcat(habitats{i}, '_converted.mat'));
        
        if exist(modelWorkspace,'file')
            fprintf('\nWorkspace exists for %s, loading...\n', habitats{i})
            load(modelWorkspace)
        else
            error('No workspace found for %s', habitats{i})
        end
        
        %% Evaluation
        plotDir = '/stud/wendering/Masterthesis/FIGURES/RAVEN_draft';
        recMethod = strcat(habitats{i},'-HMMer10E-', thresholds{t}, '-RAVEN');
        evaluationWorkspace = fullfile(RAVENModelDir,...
            strcat(habitats{i}, '_evaluation.mat'));
        evaluateModels(models, recMethod, plotDir, evaluationWorkspace, blackList)
        
    end
    
    
end
