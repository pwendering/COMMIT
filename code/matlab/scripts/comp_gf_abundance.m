%% correlate gap-filling orders and biomasses to abundances in experiments
options

modelDir = outDir;

% load medium that has been used for gap filling
load(mediumFile)

for i=1:numel(experiments)
    disp(experiments{i})
    
    % load gap-filled models
    load(fullfile(modelDir, experiments{i}));
    n = numel(GF);
    
    % load the abundances for each OTU 
    otuFile = fullfile(otuDir, habitat, '/',experiments{i},'/otutab.txt');
    otuTab = readAbundancesFromFile(otuFile);
    
    % sort by the order of models in GF
    model_ids = cellfun(@(x)strtok(x.id, '_'), GF,...
        'UniformOutput', false);
    otuTab = sortrows(otuTab, 'Row');
    [~, ia] = sortrows(otuTab, 'abundances', 'descend');
    
    % ~~~ order vs abundance ~~~ %
    
    [rho, pval] = corr(ia, gf_order', 'type', 'Spearman');
    
    fprintf('Spearman rank correlation abundance, order: %.2f (p=%.2f)\n', rho, pval)
    
    % ~~~ number of added reactions vs abundance ~~~ %
    
    added = cell(n, 1);
    for j=1:n
        % find added reactions
        added{j} = sum(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
    end
    
    [rho, pval] = corr(ia, cell2mat(added), 'type', 'Spearman');
    
    fprintf('Spearman rank correlation abundance, #added: %.2f (p=%.2f)\n', rho, pval)
    
    % ~~~ number of exported metabolites vs abundance ~~~ %
    
    exported = cell(n, 1);
    for j=1:n
        % find exported metabolites sink reactions
        exp = regexp(GF{j}.rxns(contains(GF{j}.rxns, 'sink_')),...
            'MNXM\d+\[e]', 'match');
        exp = [exp{:}]';
        exported{j} = numel(setdiff(exp, medium));
    end
    
    [rho, pval] = corr(ia, cell2mat(exported), 'type', 'Spearman');
    
    fprintf('Spearman rank correlation abundance, #exported: %.2f (p=%.2f)\n', rho, pval)
    
    % ~~~ number of exported metabolites vs abundance ~~~ %
    
    imported = cell(n, 1);
    for j=1:n
        % find imported metabolites from exchange reactions
                imp = regexp(GF{j}.rxns(contains(GF{j}.rxns, 'EX_')),...
            'MNXM\d+\[e]', 'match');
        imp = [imp{:}]';
        imported{j} = numel(setdiff(imp, medium));
    end
    
    [rho, pval] = corr(ia, cell2mat(imported), 'type', 'Spearman');
    
    fprintf('Spearman rank correlation abundance, #imported: %.2f (p=%.2f)\n', rho, pval)
    
end
