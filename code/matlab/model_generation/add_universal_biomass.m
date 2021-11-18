% Add the universal prokaryotic biomass reaction from CarveMe models to all
% models, dependent on Gram stain information (Machado et al., 2018, Nucleic Acid Res.;
% Xavier et al., 2017, Metab. Eng.)
load('data/gap-filling/biomass-reaction/universal-biomass-reaction.mat')
% biomass_rxn_uni = biomass_rxn;
% load('data/gap-filling/biomass-reaction/gramneg-biomass-reaction.mat')
% biomass_rxn_neg = biomass_rxn;
% load('data/gap-filling/biomass-reaction/grampos-biomass-reaction.mat')
% biomass_rxn_pos = biomass_rxn;
% clear biomass_rxn
biomass_name = 'universal_Biomass_reaction';
% gram_stain_info = readtable('data/gap-filling/biomass-reaction/gram-stain-wiki.tsv',...
%     'FileType','text','ReadVariableNames',false);
habitats = {'Soil', 'Leaf', 'Root'};

%% KBase models
modelDir = 'data/models/kbase';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        % choose appropriate biomass reaction dependent on gram stain
%         gr_stain = char(gram_stain_info.(3)(ismember(gram_stain_info.(1),strtok(model.id,'_'))));
%         switch gr_stain
%             case 'grampos'
%                 biomass_rxn = biomass_rxn_pos;
%             case 'gramneg'
%                 biomass_rxn = biomass_rxn_neg;
%             case ''
%                 biomass_rxn = biomass_rxn_uni;
%         end
        
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
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_biomass.mat'));
    save(workspace, 'models')
    clear models
end
        
%% RAVEN models
modelDir = 'data/models/raven/HMMer10E-50';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, strcat(habitats{i}, '_models_COBRA_GPR.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        % choose appropriate biomass reaction dependent on gram stain
%         gr_stain = char(gram_stain_info.(3)(ismember(gram_stain_info.(1),strtok(model.id,'_'))));
%         switch gr_stain
%             case 'grampos'
%                 biomass_rxn = biomass_rxn_pos;
%             case 'gramneg'
%                 biomass_rxn = biomass_rxn_neg;
%             case ''
%                 biomass_rxn = biomass_rxn_uni;
%         end
        
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
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        model = rmfield(model, 'grRules');
        models{j} = model;
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_biomass.mat'));
    save(workspace, 'models')
    clear models
end

%% CarveMe
modelDir = 'data/models/carveme';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        % choose appropriate biomass reaction dependent on gram stain
%         gr_stain = char(gram_stain_info.(3)(ismember(gram_stain_info.(1),strtok(model.id,'_'))));
%         switch gr_stain
%             case 'grampos'
%                 biomass_rxn = biomass_rxn_pos;
%             case 'gramneg'
%                 biomass_rxn = biomass_rxn_neg;
%             case ''
%                 biomass_rxn = biomass_rxn_uni;
%         end
        
        model = addReaction(model, 'BIOMASS_Reaction',...
            'reactionName', biomass_name,...
            'reactionFormula',  biomass_rxn{:},...
            'reversible', false,...
            'lowerBound', 0,...
            'upperBound', 1000, ...
            'objectiveCoef', 1,...
            'printLevel', 0);
        % add BIOMASS to the biomass reaction, add metabolite of not
        % present in the model
        if ~ismember('BIOMASS[c]',model.mets)
            model = addMetabolite(model,'BIOMASS[c]');
        end
        model = changeRxnMets(model, 'BIOMASS[c]', 'BIOMASS[c]',...
            'BIOMASS_Reaction', 1);
        % add an exchange reaction for BIOMASS
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_biomass.mat'));
    save(workspace, 'models')
    clear models
    
end

%% AuReMe
modelDir = 'data/models/aureme';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        % choose appropriate biomass reaction dependent on gram stain
%         gr_stain = char(gram_stain_info.(3)(ismember(gram_stain_info.(1),strtok(model.id,'_'))));
%         switch gr_stain
%             case 'grampos'
%                 biomass_rxn = biomass_rxn_pos;
%             case 'gramneg'
%                 biomass_rxn = biomass_rxn_neg;
%             case ''
%                 biomass_rxn = biomass_rxn_uni;
%         end
        
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
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        
        models{j} = model;
    end
    workspace = fullfile(modelDir, habitats{i}, strcat(habitats{i}, '_models_biomass.mat'));
    save(workspace, 'models')
    clear models
end

%% Merged models with CarveMe models
modelDir =  'data/models/consensus';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models.mat'));
    load(workspace)
    for j=1:numel(merged_models)
        model = merged_models{j};
        % choose appropriate biomass reaction dependent on gram stain
%         gr_stain = char(gram_stain_info.(3)(ismember(gram_stain_info.(1),strtok(model.id,'_'))));
%         switch gr_stain
%             case 'grampos'
%                 biomass_rxn = biomass_rxn_pos;
%             case 'gramneg'
%                 biomass_rxn = biomass_rxn_neg;
%             case ''
%                 biomass_rxn = biomass_rxn_uni;
%         end
        
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
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        
        merged_models{j} = model;
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models_biomass.mat'));
    save(workspace, 'merged_models')
    clear merged_models
end

%% Merged models without CarveMe models

modelDir = 'data/models/consensus/';
for i=1:numel(habitats)
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models_noCarveMe.mat'));
    load(workspace)
    for j=1:numel(merged_models)
        model = merged_models{j};
        % choose appropriate biomass reaction dependent on gram stain
%         gr_stain = char(gram_stain_info.(3)(ismember(gram_stain_info.(1),strtok(model.id,'_'))));
%         switch gr_stain
%             case 'grampos'
%                 biomass_rxn = biomass_rxn_pos;
%             case 'gramneg'
%                 biomass_rxn = biomass_rxn_neg;
%             case ''
%                 biomass_rxn = biomass_rxn_uni;
%         end
        
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
        [~,model]=evalc('addSinkReactions(model, ''BIOMASS[c]'', 0, 1000);');
        
        merged_models{j} = model;
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_consensus_models_noCarveMe_biomass.mat'));
    save(workspace, 'merged_models')
    clear merged_models
end
