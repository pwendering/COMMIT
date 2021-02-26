% compare distance measures of merged models to sequence similarity
options

tablesDir = fullfile(topDir, 'data/tables');
figOutDir = fullfile(topDir, 'figures/Comparison-of-methods');

% Load distance matrix from newick tree
phyloFile = fullfile(topDir, 'data/genomes/Phylogeny/nw_distance_AtSPHERE.txt');
phylo_dist = readtable(phyloFile,...
    'ReadVariableNames', true,...
    'ReadRowNames', true);

ecTranslationTable = readtable(fullfile(tablesDir, 'corrected-EC-numbers.csv'),...
    'ReadVariableNames', false);

coFactorsTab = readtable(fullfile(tablesDir, 'cofactors_from_KEGG.csv'),...
    'Delimiter', '\t', 'ReadVariableNames', false);
coFactors = coFactorsTab.Var5; % MNXref

habitats = {'Soil', 'Root', 'Leaf'};
labels = {'SVD', 'JD_r', 'JD_m', 'n_DE', 'JD_DE', 'JD_EC', 'c_EC', 'c_CF'};

parpool(ncpus);
for i=1:numel(habitats)
    %% Load workspaces
    disp(habitats{i})
    % Consensus models including CarveMe reconstructions
    models = load(fullfile(topDir, 'data/models/consensus/',...
        [habitats{i}, '_consensus_models.mat']), 'merged_models');
    models = models.merged_models;
    %     merged_noCarveMe = load(fullfile(topDir, 'data/models/consensus/',...
    %         [habitats{i}, '_consensus_models_noCarveMe.mat']), 'merged_models');
    %     merged_noCarveMe = merged_noCarveMe.merged_models;
    
    
    %% calculate distances
    % 1. SVD distance
    SVD_matrix = svdDistance(models);
    % 2. Reaction Jaccard distance
    rxnJaccard_matrix = rxnJaccardDist(models);
    % 3. metabolite Jaccard distance
    metJaccard_matrix = metJaccardDist(models);
    % 4. number of dead-ends distance
    deadEndsNumber_matrix = deadEndsNumberDist(models);
    % 5. dead-end metabolite Jaccard distance
    deadEndJaccard_matrix = deadEndsJaccardDist(models);
    % 6. E.C. number Jaccard distance (all four levels)
    ecLevel = 4;
    ecJaccard_matrix = ecJaccardDist(models, ecLevel);
    % 7. abundance of E.C. numbers distance
    ecAbundanceDist_matrix = ecAbundanceDist(models, ecTranslationTable.Var2);
    % 8. co-factor usage distance
    coFactorDistance_matrix = coFactorDist(models, coFactors);
    
    M = {SVD_matrix, rxnJaccard_matrix, metJaccard_matrix, deadEndsNumber_matrix,...
        deadEndJaccard_matrix, ecJaccard_matrix, ecAbundanceDist_matrix,...
        coFactorDistance_matrix};
    
    %% Compare distance matrices to phylogeny
    
    % re-order phylogeny matrix a cording to models array
    model_ids = cellfun(@(x)x.id, models, 'UniformOutput', false);
    
    tax_names = cell(numel(models),1);
    tax_names(:) = {''};
    for j=1:numel(models)
        name = models{j}.description;
        split = strsplit(name ,'_');
        count = 1;
        name = strcat(split{1},'_', num2str(count));
        while ismember(name, tax_names)
            count = count + 1;
            name = regexpi(name, '[a-z]*', 'match');
            name = strcat(name{:}, '_',  num2str(count));
        end
        tax_names{j} = name;
    end
    clear split name count


    row_names = phylo_dist.Properties.RowNames;
    phylo_dist_habitat = phylo_dist(contains(row_names, habitats{i}), contains(row_names, habitats{i}));
    row_names = phylo_dist_habitat.Properties.RowNames;
    re_order = cellfun(@(x)find(strcmp(x, row_names)), model_ids);
    phylo_dist_habitat = table2array(phylo_dist_habitat(re_order, re_order)); clear re_order row_names
    
    mantel_corr_dist_phyl = zeros(size(M));
    p_emp = zeros(size(M));
    
    
    for j=1:numel(M)
        matrix = M{j};
        % calculate pearson correlation 
        mantel_corr_dist_phyl(j) = mantelCoefficient(phylo_dist_habitat, matrix);
        
        % assess significance
        c = 0;
        iter = 10000;
        for k=1:iter
            [~,ia] = sort(rand(size(matrix)));
            R = matrix(ia);
            r = mantelCoefficient(phylo_dist_habitat, R);
            c = c + (r > mantel_corr_dist_phyl(j));
        end
        p_emp(j) = (c + 1) / (iter + 1);
        
        D = array2table(M{j},...
            'VariableNames', model_ids,...
            'RowNames', tax_names);
        
        writetable(D, fullfile(figOutDir, [labels{j},'-', habitats{i}, '.txt']),...
            'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t')
    end
    
    
    writetable(cell2table([labels', num2cell(mantel_corr_dist_phyl'), num2cell(p_emp')],...
        'VariableNames', {'Measure', 'Mantel_r', 'p_value'}),...
        fullfile(figOutDir, ['dist-phylo-', habitats{i}, '.txt']),...
        'WriteVariableNames', true,...
        'Delimiter', '\t')
    
    
end
