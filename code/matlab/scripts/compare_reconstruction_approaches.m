% Compare the reconstruction approaches with each other and with the
% consensus models (with and wihtout CarveMe)

% load options to get top directory
options

tablesDir = 'data/tables';

% specify workspace file name where figure data will be saved
figDataDir = 'figures/Comparison-of-methods';

ecTranslationTable = readtable(fullfile(tablesDir, 'corrected-EC-numbers.csv'),...
    'ReadVariableNames', false);

coFactorsTab = readtable(fullfile(tablesDir, 'cofactors_from_KEGG.csv'),...
    'Delimiter', '\t', 'ReadVariableNames', false);
coFactors = coFactorsTab.Var5; % MNXref

n_models = 432;
habitats = {'Soil', 'Leaf', 'Root'};
OTU_dist_matrices = cell(n_models, 1);
geneJaccard_matrices = cell(n_models, 1);
c = 0;

parpool(2);
for i=1:numel(habitats)
    %% Load workspaces
    disp(habitats{i})
    % KBase models
    KB_models = load(fullfile('data/models/kbase', habitats{i},...
        [habitats{i}, '_models_no_medium_no_biomass.mat']), 'models');
    KB_models = KB_models.models;
    % CarveMe models
    CM_models = load(fullfile('data/models/carveme', habitats{i},...
        [habitats{i}, '_models_no_medium_no_biomass.mat']), 'models');
    CM_models = CM_models.models;
    % RAVEN 2.0 models
    RV_models = load(fullfile('data/models/raven','HMMer10E-50',...
        [habitats{i}, '_models_no_medium_no_biomass.mat']), 'models');
    RV_models = RV_models.models;
    % AuReMe models
    AU_models = load(fullfile('data/models/aureme',...
        habitats{i},[habitats{i}, '_models_no_medium_no_biomass.mat']), 'models');
    AU_models = AU_models.models;
    % Consensus models including CarveMe reconstructions
    merged = load(fullfile('data/models/consensus',...
        [habitats{i}, '_consensus_models.mat']), 'merged_models');
    merged = merged.merged_models;
    %     merged_noCarveMe = load(fullfile('data/models/consensus',...
    %         [habitats{i}, '_consensus_models_noCarveMe.mat']), 'merged_models');
    %     merged_noCarveMe = merged_noCarveMe.merged_models;
    for j=1:numel(KB_models)
        
        c = c + 1;
        
        models = {KB_models{i}, CM_models{i}, RV_models{i}, AU_models{i},...
            merged{i}};
        
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
        
        OTU_dist_matrices{c} = STATIS(M);
        
        geneJD_mt = zeros(numel(models));
        
        for k=1:numel(models)-1
            for l=k+1:numel(models)
                geneJD_mt(k,l) = 1 - ...
                    numel(intersect(models{k}.genes, models{l}.genes)) / ...
                    numel(union(models{k}.genes, models{l}.genes));
                geneJD_mt(l,k) = geneJD_mt(k,l);
            end
        end
        
        geneJaccard_matrices{c} = geneJD_mt;
        
    end
    
end; clear *models merged ecTranslationTable coFactor* *matrix


%% For each OTU, compare the four individual models using distance measures and combine them into a single matrix
method_dist_mat = STATIS(OTU_dist_matrices);
gene_dist_mat = STATIS(geneJaccard_matrices);
labels = {'KBase', 'CarveMe', 'RAVEN', 'AuReMe', 'Consensus'};
save(fullfile(figDataDir, 'fig-data'), 'gene_dist_mat', 'method_dist_mat')

load(fullfile(figDataDir, 'fig-data'))
dist_mt = method_dist_mat;
for i=2:size(dist_mt,2)
    for j=1:i-1
        dist_mt(i,j) = gene_dist_mat(i,j);
    end
end

writetable(array2table(dist_mt,...
    'VariableNames', labels, 'RowNames', labels),...
    fullfile(figDataDir, 'dist_mt.txt'))



