% compare numbers of reactions, metabolites and genes in the merged models
% compared to the single draft models

% load options to get top directory
options

habitats = {'Soil', 'Leaf', 'Root'};
methods = {'CarveMe', 'KBase', 'AuReMe', 'RAVEN', 'consensus'};
features = {'rxns', 'mets', 'genes'};

modelDir = fullfile('data','models');
figOutDir = fullfile('figures','model-features');

for i=1:numel(habitats)
    clear sums*
    disp(habitats{i})
    disp('-----------------')
    % initialize array for numbers of reactions with gene association
    % across the reconstruction approaches
    perc_rxns_with_gpr = zeros(500,numel(methods));
    
    for j=1:numel(methods)
        fprintf('\t> %s\n', methods{j})
        switch methods{j}
            case 'RAVEN'
                workspace = fullfile(modelDir, lower(methods{j}),...
                    'HMMer10E-50', strcat(habitats{i}, '_models_COBRA_GPR'));
            case 'consensus'
                workspace = fullfile(modelDir, lower(methods{j}), strcat(habitats{i},...
                    '_consensus_models.mat'));
            otherwise
                workspace = fullfile(modelDir, lower(methods{j}),...
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
        
        % Numbers of features
        for k=1:numel(features)
            eval([methods{j}, '_', features{k},...
                ' = cell2mat(cellfun(@(x)numel(x.', features{k}, '), models, ''UniformOutput'', false));']);
            if ~isequal(methods{j}, 'consensus')
                eval(['sums_', features{k}, ' = sums_', features{k}, ' + ', methods{j}, '_', features{k},';']);
            end
        end
        
        % Rreactions with gene association
        % KBase models contain a genes "{''}", which appears in multiple GPR
        % rules
        models = cellfun(@(M)removeGenesFromModel(M,{''}),models,'un',0);
        for k=1:numel(models)
           perc_rxns_with_gpr(k,j) = sum(~cellfun(@isempty,models{k}.rules))/numel(models{k}.rxns);
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
    
    perc_rxns_with_gpr(all(perc_rxns_with_gpr==0,2),:) = [];
    writetable(array2table(perc_rxns_with_gpr,'VariableNames',methods),...
        fullfile(figOutDir, [habitats{i} '_rxns_gpr_percent.txt']),...
        'Delimiter', '\t')
    
    disp('')
end
