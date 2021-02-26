% compare numbers of reactions, metabolites and genes in the merged models
% compared to the single draft models

% load options to get top directory
options
clearvars -except topDir

habitats = {'Soil', 'Leaf', 'Root'};
methods = {'CarveMe', 'KBase', 'AuReMe', 'RAVEN', 'consensus'};
features = {'rxns', 'mets', 'genes'};

dataDir = fullfile(topDir, 'data');
figOutDir = fullfile(topDir, 'figures/model-features');

for i=1:numel(habitats)
    clear sums*
    disp(habitats{i})
    disp('-----------------')
    for j=1:numel(methods)
        fprintf('\t> %s\n', methods{j})
        switch methods{j}
            case 'RAVEN'
                workspace = fullfile(dataDir, strcat('models_', methods{j}),...
                    'HMMer10E-50', strcat(habitats{i}, '_models_COBRA_GPR'));
            case 'consensus'
                workspace = fullfile(dataDir, 'Consensus_models', strcat(habitats{i},...
                    '_consensus_models.mat'));
            otherwise
                workspace = fullfile(dataDir, strcat('models_', methods{j}),...
                    habitats{i}, strcat(habitats{i}, '_models_metFormulas'));
        end
        
        load(workspace);
        
        if ~exist('models', 'var')
            models = merged_models; clear merged_models
        end
        
        models = reshape(models, numel(models), 1);
        
        if ~exist('sums_rxns', 'var')
            sums_rxns = zeros(size(models));
            sums_mets = zeros(size(models));
            sums_genes = zeros(size(models));
        end
        
        for k=1:numel(features)
            eval([methods{j}, '_', features{k},...
                ' = cell2mat(cellfun(@(x)numel(x.', features{k}, '), models, ''UniformOutput'', false));']);
            if ~isequal(methods{j}, 'consensus')
                eval(['sums_', features{k}, ' = sums_', features{k}, ' + ', methods{j}, '_', features{k},';']);
            end
        end
        
        clear models workspace
    end
        
    % write results
    disp('==> writing results to file')
    for j=1:numel(features)
        filename = fullfile(figOutDir, [features{j}, '-', habitats{i}, '.txt']);
        eval(['outArray = [', strjoin(strcat([methods, {'sums'}],'_', features{j}),', '), '];']);
        writetable(array2table(outArray,...
            'VariableNames', strcat([methods, {'sums'}],'_', features{j})),...
            filename, 'Delimiter', '\t',...
            'WriteVariableNames', true, 'WriteRowNames', false);
    end
    disp('')
end
clear topDir