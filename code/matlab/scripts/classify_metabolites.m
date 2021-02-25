
warning('off')

% select habitat
habitat = 'Soil';

% medium that has been used for gap filling
mediumFile = '/stud/wendering/Masterthesis/DATA/media/minimal-medium.mat';
load(mediumFile)

% gap-filling database which contains permeability information
load('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat')

% select experiments to compare
experiments = {'Schlaeppi', 'Bulgarelli'};
subFolder = 'all';

% taxonomic classification of At-SPHERE OTUs
% taxonomyFile = '/stud/wendering/Masterthesis/DATA/genomes/At-SPHERE-phyla.txt';
% taxonomyFile = '/stud/wendering/Masterthesis/DATA/genomes/At-SPHERE-classes.txt';
taxonomyFile = '/stud/wendering/Masterthesis/DATA/genomes/At-SPHERE-families.txt';

% directories for metabolic models to compare
iterativeGfDir = fullfile('/stud/wendering/Masterthesis/DATA/Gap-filling/iterative',...
    habitat, subFolder);

% load DAG for ChEBI ontology terms and associated identifiers for every
% node
load('/net/calc1/srv/wendering/chebi/ontologyGraph.mat', 'ontologyMat')

% tab-separated ontology file
ontologyFile = '/net/calc1/srv/wendering/chebi/chebi_ontology.csv';
ontologyTab = readtable(ontologyFile, 'ReadVariableNames', true,...
    'Delimiter', '\t');
names = ontologyTab.name;
terms = strtrim(cellstr(num2str(ontologyTab.ID)));
clear ontologyTab


%% Contingency table for exported and imported metabolites
for i=1:numel(experiments)
    
    disp(experiments{i})
    
    % read gap filled consensus models
    load(fullfile(iterativeGfDir, experiments{i}), 'GF')
    
    model_ids = cellfun(@(x)strtok(x.id, '_'), GF,...
        'UniformOutput', false);
    
    % extract exported metabolites from sink reactions that have been added
    % suring iterative gap filling
    exported_per_model = cell(numel(GF), 1);
    imported_per_model = cell(numel(GF), 1);
    
    for j=1:numel(GF)
        
        % find exported metabolites sink reactions
        exported = regexp(GF{j}.rxns(contains(GF{j}.rxns, 'sink_')),...
            'MNXM\d+\[e]', 'match');
        exported_per_model{j} = [exported{:}]'; clear exported
        
        % find imported metabolites from exchange reactions
        imported = regexp(GF{j}.rxns(contains(GF{j}.rxns, 'EX_')),...
            'MNXM\d+\[e]', 'match');
        imported = [imported{:}]';
        imported_per_model{j} = setdiff(imported, medium);
        
    end; clear GF 
    
    % get permeability for all metabolites and determine the ones that are not
    % predicted to be permeable
    
    % imported
    all_imported = vertcat(imported_per_model{:});
    all_imported = unique(all_imported);
        
    p_imported = ...
        cellfun(@(x)dbModel_MNXref_balanced.metPermeability(strcmp(x,...
        dbModel_MNXref_balanced.mets)), all_imported, 'un', 0);
    p_imported(cellfun('isempty', p_imported)) = {false};
    p_imported = cell2mat(p_imported);
    mets_not_p = all_imported(~p_imported);
    imported_per_model = cellfun(@(x)setdiff(x, mets_not_p), imported_per_model,...
        'UniformOutput', 0);
    
    % exported
    all_exported = vertcat(exported_per_model{:});
    all_exported = unique(all_exported);

    p_exported = ...
        cellfun(@(x)dbModel_MNXref_balanced.metPermeability(strcmp(x,...
        dbModel_MNXref_balanced.mets)), all_exported, 'un', 0);
    p_exported(cellfun('isempty', p_exported)) = {false};
    p_exported = cell2mat(p_exported);
    mets_not_p = all_exported(~p_exported);
    exported_per_model = cellfun(@(x)setdiff(x, mets_not_p), exported_per_model,...
        'UniformOutput', 0);
    
    % merge results per taxonomic class
    taxonomyTab = readtable(taxonomyFile, 'readVariableNames', true, 'Delimiter', '\t');
    taxonomyTab = taxonomyTab(ismember(taxonomyTab.isolate_ID,model_ids),:);
    tax_classes = unique(taxonomyTab.(2));
    
    exported_per_class = cell(numel(tax_classes), 1);
    imported_per_class = cell(numel(tax_classes), 1);
    
    abundances_per_class_exp = zeros(numel(tax_classes), numel(terms));
    abundances_per_class_imp = zeros(numel(tax_classes), numel(terms));
    
    fprintf('\t> counting abundances...\n')
    exclude = 3;
    %% Enrichment analysis
    for j=1:numel(tax_classes)
        class_id_matching = ismember(taxonomyTab.(2), tax_classes{j});
        exported_per_class{j} = unique(vertcat(exported_per_model{class_id_matching}));
        imported_per_class{j} = unique(vertcat(imported_per_model{class_id_matching}));
        
        % translate metabolite IDs to ChEBI namespace
        metList_exp = translateIDs(strtok(exported_per_class{j}, '['), 'met',...
            [], 'MNXref', 'ChEBI', false);
        metList_imp = translateIDs(strtok(imported_per_class{j}, '['), 'met',...
            [], 'MNXref', 'ChEBI', false);
        
        % include only the entries that start with a number i.e. ChEBI ID
        metList_exp = regexp(metList_exp, '^\d.*', 'match');
        metList_exp = [metList_exp{:}]';
        metList_imp = regexp(metList_imp, '^\d.*', 'match');
        metList_imp = [metList_imp{:}]';
        
        % for every metabolite choose the first ChEBI identifier that is
        % contained in the list of all identifiers/ontology terms
        for k=1:numel(metList_exp)
            tmp_ids = regexp(metList_exp{k}, '\d+', 'match');
            tmp_idx = find(ismember(tmp_ids, terms));
            metList_exp{k} = tmp_ids{tmp_idx(1)};
        end
        
        for k=1:numel(metList_imp)
            tmp_ids = regexp(metList_imp{k}, '\d+', 'match');
            tmp_idx = find(ismember(tmp_ids, terms));
            metList_imp{k} = tmp_ids{tmp_idx(1)};
        end
        
        exported_per_class{j} = metList_exp;
        imported_per_class{j} = metList_imp;
        
        % determine of all terms in the given test set
        abundances_per_class_exp(j,:) = getAbundances(ontologyMat, metList_exp, terms, exclude);
        abundances_per_class_imp(j,:) = getAbundances(ontologyMat, metList_imp, terms, exclude);
        
    end
    
    
    % remove empty columns and terms from the metabolite sets...
    idx_nz_exp = any(abundances_per_class_exp,1);
    idx_nz_exp(ismember(terms, vertcat(exported_per_class{:}))) = 0;
    idx_nz_imp = any(abundances_per_class_imp,1);
    idx_nz_imp(ismember(terms, vertcat(imported_per_class{:}))) = 0;
    
    % ...and cut the reverted DAG at a specific level
    level = 3;
    G_undirected = ontologyMat + ontologyMat';
    
    idx_nz_exp(graphtraverse(G_undirected,...
        1, 'Method', 'BFS', 'Depth', level)) = 0;
    idx_nz_imp(graphtraverse(G_undirected,...
        1, 'Method', 'BFS', 'Depth', level)) = 0;
    clear G_undirected
    
    abundances_per_class_exp = abundances_per_class_exp(:, idx_nz_exp);
    abundances_per_class_imp = abundances_per_class_imp(:, idx_nz_imp);
    
    names_exp = names(idx_nz_exp);
    names_imp = names(idx_nz_imp);
    
    n_met_classes_exp = size(abundances_per_class_exp,2);
    n_met_classes_imp = size(abundances_per_class_imp,2);
    
    
    %% Testing
    
    fprintf('\t> one-sided hypergeometric test for every term...\n')
    alpha = 0.05;
    % for each metabolic class (MC):
    %       compare the abundance in each taxonomic class (TC) to the
    %       abundance over all TCs
    % Two-sided hypergeometric test
    % e.g.
    %         MC1 | not MC1
    %        --------------
    % TC1      8  |  1472   | 1480
    % TC2-3    16 |  2085   | 2101
    %        ----------------------
    %          24 |  3557   | 3581
    
    % MATLAB implemetation for the two-sided test close to minimum
    % likelyhood p-value
    %     p_two_sided = zeros(numel(tax_classes), n_met_classes_exp);
    p_one_sided_upper_exp = zeros(numel(tax_classes), n_met_classes_exp);
    %     p_one_sided_lower = zeros(numel(tax_classes), n_met_classes_exp);
    
    %     X = zeros(2);
    M = sum(sum(abundances_per_class_exp,1));
    
    for tax=1:numel(tax_classes)
        for term=1:n_met_classes_exp
            
            % contigency table for current term
            %             X(1,1) = abundances_per_class_exp(tax, term);
            %             X(2,1) = sum(abundances_per_class_exp(setdiff(1:numel(tax_classes), tax), term));
            %             X(1,2) = sum(abundances_per_class_exp(tax, setdiff(1:size(abundances_per_class_exp, 2), term)));
            %             X(2,2) = sum(sum(abundances_per_class_exp(setdiff(1:numel(tax_classes), tax), :), 1)) - X(2,1);
            %             [~, p_two_sided(tax,term)] = fishertest(X);
            
            x = abundances_per_class_exp(tax, term);
            K = sum(abundances_per_class_exp(:,term));
            N = sum(abundances_per_class_exp(tax, :));
            
            p_one_sided_upper_exp(tax,term) = hygecdf(x,M,K,N, 'upper');
            %             p_one_sided_lower(tax,term) = hygecdf(x,M,K,N);
        end
    end
    
    
    p_one_sided_upper_imp = zeros(numel(tax_classes), n_met_classes_imp);
    %     p_one_sided_lower = zeros(numel(tax_classes), n_met_classes_exp);
    
    %     X = zeros(2);
    M = sum(sum(abundances_per_class_imp,1));
    
    for tax=1:numel(tax_classes)
        for term=1:n_met_classes_imp
            
            x = abundances_per_class_imp(tax, term);
            K = sum(abundances_per_class_imp(:,term));
            N = sum(abundances_per_class_imp(tax, :));
            
            p_one_sided_upper_imp(tax,term) = hygecdf(x,M,K,N, 'upper');
        end
    end
    
    
    %% multiple testing correction
    
    % FDR
    
    % mafdr function in MATLAB
    %     h_two_sided = mafdr(reshape(p_two_sided, numel(p_two_sided), 1));
    %     [~, ranks] = sort(reshape(p_two_sided, numel(p_two_sided), 1), 'ascend');
    %
    %     p_fdr = ranks/M *alpha;
    %     h = reshape(p_two_sided, numel(p_two_sided), 1) <= p_fdr;
    %     h = reshape(h, numel(tax_classes), n_met_classes);
    
    %     mafdr(reshape(p_one_sided_upper, numel(p_one_sided_upper), 1))
    
    [p_sorted, idx] = sort(reshape(p_one_sided_upper_exp, numel(p_one_sided_upper_exp), 1), 'ascend');
    
    p_fdr = (1:numel(idx))/numel(idx) *alpha;
    p_sorted = p_sorted <= p_fdr';
    h_one_sided_upper_exp = reshape(p_sorted(idx), numel(tax_classes), n_met_classes_exp);
    
    [p_sorted, idx] = sort(reshape(p_one_sided_upper_imp, numel(p_one_sided_upper_imp), 1), 'ascend');
    
    p_fdr = (1:numel(idx))/numel(idx) *alpha;
    p_sorted = p_sorted <= p_fdr';
    h_one_sided_upper_imp = reshape(p_sorted(idx), numel(tax_classes), n_met_classes_imp);
    
    fprintf('\t> finding enriched terms per clade...\n')
    
    %% get enriched terms per clade
    enriched_per_clade_exp = cell(numel(tax_classes), 1);
    for j=1:numel(tax_classes)
        enriched_per_clade_exp{j} = names_exp(h_one_sided_upper_exp(j,:));
    end
    
    enriched_per_clade_imp = cell(numel(tax_classes), 1);
    for j=1:numel(tax_classes)
        enriched_per_clade_imp{j} = names_exp(h_one_sided_upper_imp(j,:));
    end
    
    %% write results to file
    fprintf('\t> writing results...\n')
    
    fid = fopen(fullfile('/stud/wendering/Masterthesis/DATA/Gap-filling/iterative/',...
        habitat, subFolder, [experiments{i}, '-evaluation.txt']), 'w');
    
    fprintf(fid,...
        'Enriched classes among exchanged metabolites for phyla observed in %s dataset\n\n',...
        experiments{i});
    fprintf(fid,...
        ['The third neighborhood of metabolite terms was removed from the universal set.\n',...
        'The second neighborhood of the root was removed from the universal set.\n\n']);
    fprintf(fid,...
        'For each metabolite ontology term, a one-sided Fisher test for the upper tail was performed.\n\n');
    fprintf(fid,...
        'All p-values were corrected using the Benjamini-Hochberg procedure (initial alpha = %.2f).\n\n',...
        alpha);
    
%     m_exchange = zeros(numel(tax_classes));
    
    for j=1:numel(tax_classes)
        n_current_tax = sum(ismember(taxonomyTab.(2), tax_classes{j}));
        fprintf(fid, '%s (%d)\n\n', tax_classes{j}, n_current_tax);
        fprintf(fid, 'Exported (%d):\n\n', numel(enriched_per_clade_exp{j}));
        fprintf(fid, '%s\n', enriched_per_clade_exp{j}{:});
        fprintf(fid, '\n');
        fprintf(fid, 'Imported(%d):\n\n', numel(enriched_per_clade_imp{j}));
        fprintf(fid, '%s\n', enriched_per_clade_imp{j}{:});
        fprintf(fid, '\n');
        
        % pairwise intersection of export and import between classes
%         for k=1:numel(tax_classes)
%             if j~=k
%                 m_exchange(j,k) = numel(intersect(enriched_per_clade_exp{j},...
%                     enriched_per_clade_imp{k}));
%             end
%         end
    end
%     fprintf(fid, 'exchange matrix: pairwise intersection of import and export\n');
%     pattern = [strjoin(repmat({'%d'}, 1, numel(tax_classes)), '\\t'), '\n'];
%     for j=1:size(m_exchange, 2)
%         fprintf(fid, pattern, m_exchange(j,:));
%     end
    fclose(fid);
    
    
end
