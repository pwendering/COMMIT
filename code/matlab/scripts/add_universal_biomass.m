% Add the universal prokaryotic biomass reaction from CarveMe models to all
% models (Machado et al., 2018, Nucleic Acid Res.; Xavier et al., 2017, Metab. Eng.)
load('/stud/wendering/Masterthesis/DATA/Gap-filling/universal-biomass-reaction.mat')
biomass_name = 'universal_Biomass_reaction'; 
habitats = {'Soil', 'Leaf', 'Root'};

%% KBase models
modelDir = '/stud/wendering/Masterthesis/DATA/models_KBase';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        model = addReaction(model, 'BIOMASS_Reaction',...
            'reactionName', biomass_name,...
            'reactionFormula',  biomass_rxn{:},...
            'reversible', false,...
            'lowerBound', 0,...
            'upperBound', 1000, ...
            'objectiveCoef', 1,...
            'printLevel', 0);
        % add BIOMASS as a product to the biomass reaction
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        % add an exchange reaction for BIOMASS
        model = addSinkReactions(model, 'BIOMASS[c]', 0, 1000);
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_biomass.mat'));
    save(workspace, 'models')
    clear models
end
        
%% RAVEN models
thresholds = {'50'};%, '30'};

for t = thresholds
    modelDir = fullfile('/stud/wendering/Masterthesis/DATA/models_RAVEN', ...
        strcat('HMMer10E-', t{:}));
    for i=1:numel(habitats)
        workspace = fullfile(modelDir, strcat(habitats{i}, '_models_COBRA_GPR.mat'));
        load(workspace)
        for j=1:numel(models)
            model = models{j};
            model = addReaction(model, 'BIOMASS_Reaction',...
                'reactionName', biomass_name,...
                'reactionFormula',  biomass_rxn{:},...
                'reversible', false,...
                'lowerBound', 0,...
                'upperBound', 1000, ...
                'objectiveCoef', 1,...
                'printLevel', 0);
            model.c(end) = 1;
            % add BIOMASS as a metabolite
            model = addMetabolite(model, 'BIOMASS[c]', 'metName', 'Biomass');
            % add BIOMASS to the biomass reaction
            model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
                'BIOMASS_Reaction', 1);
            % add an exchange reaction for BIOMASS
            model = addSinkReactions(model, 'BIOMASS[c]', 0, 1000);
            model = rmfield(model, 'grRules');
            models{j} = model;
        end
        workspace = fullfile(modelDir, strcat(habitats{i}, '_biomass.mat'));
        save(workspace, 'models')
        clear models
    end
    
end

%% CarveMe
modelDir = '/stud/wendering/Masterthesis/DATA/models_CarveMe';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        model = addReaction(model, 'BIOMASS_Reaction',...
                'reactionName', biomass_name,...
                'reactionFormula',  biomass_rxn{:},...
                'reversible', false,...
                'lowerBound', 0,...
                'upperBound', 1000, ...
                'objectiveCoef', 1,...
                'printLevel', 0);
        % add BIOMASS to the biomass reaction
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        % add an exchange reaction for BIOMASS
        model = addSinkReactions(model, 'BIOMASS[c]', 0, 1000);
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_biomass.mat'));
    save(workspace, 'models')
    clear models

end

%% AuReMe
modelDir = '/stud/wendering/Masterthesis/DATA/models_AuReMe';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        model = addReaction(model, 'BIOMASS_Reaction',...
            'reactionName', biomass_name',...
            'reactionFormula',  biomass_rxn{:},...
            'reversible', false,...
            'lowerBound', 0,...
            'upperBound', 1000, ...
            'objectiveCoef', 1,...
            'printLevel', 0);
        % add BIOMASS as a metabolite
        model = addMetabolite(model, 'BIOMASS[c]', 'metName', 'Biomass');
        % add BIOMASS as a product to the biomass reaction
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        % add an exchange reaction for BIOMASS
        model = addSinkReactions(model, 'BIOMASS[c]', 0, 1000);
        
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_biomass.mat'));
    save(workspace, 'models')
    clear models
end

%% Merged models with CarveMe models
modelDir = '/stud/wendering/Masterthesis/DATA/Consensus_models';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models.mat'));
    load(workspace)
    for j=1:numel(merged_models)
        model = merged_models{j};
        model = addReaction(model, 'BIOMASS_Reaction',...
            'reactionName', biomass_name',...
            'reactionFormula',  biomass_rxn{:},...
            'reversible', false,...
            'lowerBound', 0,...
            'upperBound', 1000, ...
            'objectiveCoef', 1,...
            'printLevel', 0);

        % add BIOMASS as a product to the biomass reaction
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        model = removeRxns(model, 'EX_BIOMASS_c');
        % add an exchange reaction for BIOMASS
        model = addSinkReactions(model, 'BIOMASS[c]', 0, 1000);
        
        merged_models{j} = model;
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models_biomass.mat'));
    save(workspace, 'merged_models')
    clear merged_models
end

%% Merged models without CarveMe models

modelDir = '/stud/wendering/Masterthesis/DATA/Consensus_models';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models_noCarveMe.mat'));
    load(workspace)
    for j=1:numel(merged_models)
        model = merged_models{j};
        model = addReaction(model, 'BIOMASS_Reaction',...
            'reactionName', biomass_name',...
            'reactionFormula',  biomass_rxn{:},...
            'reversible', false,...
            'lowerBound', 0,...
            'upperBound', 1000, ...
            'objectiveCoef', 1,...
            'printLevel', 0);

        % add BIOMASS as a product to the biomass reaction
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        model = removeRxns(model, 'EX_BIOMASS_c');
        % add an exchange reaction for BIOMASS
        model = addSinkReactions(model, 'BIOMASS[c]', 0, 1000);
        
        merged_models{j} = model;
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models_noCarveMe_biomass.mat'));
    save(workspace, 'merged_models')
    clear merged_models
end
