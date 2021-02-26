% find exchanged metabolites between the OTUs
options

experiment = 'Bulgarelli';

% medium that has been used for gap filling
load(mediumFile)

% gap-filling database which contains permeability information
load(dbFile)

% taxonomic classification of At-SPHERE OTUs
% taxonomyFile = fullfile(topDir, 'data/genomes/At-SPHERE-phyla.txt');
taxonomyFile = fullfile(topDir, 'data/genomes/At-SPHERE-classes.txt');
% taxonomyFile = fullfile(topDir, 'data/genomes/At-SPHERE-families.txt');

% load model workspace
modelWorkspace = fullfile(outDir, experiment);
load(modelWorkspace, 'GF')
n = numel(GF);

% read OTU OTU abundances from file
otutab = readAbundancesFromFile(fullfile(otuDir, habitat, experiment, 'otutab.txt'));
otutab = sortrows(otutab, 'Row');

%% extract exported metabolites from sink reactions that have been added
% during iterative gap filling
exported_per_model = cell(n, 1);
imported_per_model = cell(n, 1);

for i=1:n
    
    % find exported metabolites sink reactions
    exported = regexp(GF{i}.rxns(contains(GF{i}.rxns, 'sink_')),...
        'MNXM\d+\[e]', 'match');
    exported = [exported{:}]';
    exported_per_model{i} = setdiff(exported, medium); clear exported
    
    % find imported metabolites from exchange reactions
    imported = regexp(GF{i}.rxns(contains(GF{i}.rxns, 'EX_')),...
        'MNXM\d+\[e]', 'match');
    imported = [imported{:}]';
    imported_per_model{i} = setdiff(imported, medium); clear imported
    
end

%% exclude metabolites that appear in all models
t_unspecific = n;

% imported
all_imported = vertcat(imported_per_model{:});
all_imported = unique(all_imported);

% get permeability for all metabolites and determine the ones that are not
% predicted to be permeable
p_imported = ...
    cellfun(@(x)dbModel_MNXref_balanced.metPermeability(strcmp(x,...
    dbModel_MNXref_balanced.mets)), all_imported, 'un', 0);
p_imported(cellfun('isempty', p_imported)) = {false};
p_imported = cell2mat(p_imported);
mets_not_p = all_imported(~p_imported);
imported_per_model = cellfun(@(x)setdiff(x, mets_not_p), imported_per_model,...
    'UniformOutput', 0);
% exclude imported metabolites that occur in most or all models
unspecific_imported = cellfun(@(x)sum(contains(all_imported, x)), unique(all_imported))==t_unspecific;
unspecific_imported = all_imported(unspecific_imported);
imported_per_model = cellfun(@(x)setdiff(x, unspecific_imported), imported_per_model,...
    'UniformOutput', 0);

% exported
all_exported = vertcat(exported_per_model{:});
all_exported = unique(all_exported);
% get permeability for all metabolites and determine the ones that are not
% predicted to be permeable
p_exported = ...
    cellfun(@(x)dbModel_MNXref_balanced.metPermeability(strcmp(x,...
    dbModel_MNXref_balanced.mets)), all_exported, 'un', 0);
p_exported(cellfun('isempty', p_exported)) = {false};
p_exported = cell2mat(p_exported);
mets_not_p = all_exported(~p_exported);
exported_per_model = cellfun(@(x)setdiff(x, mets_not_p), exported_per_model,...
    'UniformOutput', 0);

% exclude exported metabolites that occur in most or all models
unspecific_exported = cellfun(@(x)sum(contains(all_exported, x)), unique(all_exported))==t_unspecific;
unspecific_exported = all_exported(unspecific_exported);
exported_per_model = cellfun(@(x)setdiff(x, unspecific_exported), exported_per_model,...
    'UniformOutput', 0);

clear mets_not_p all_exported all_imported unspecific_exported unspecific_imported ...
    dbModel_MNXref_balanced

%% classify metabolites by KEGG pathway and brite definitions per model
model_ids = cellfun(@(x)strtok(x.id, '_'), GF,...
    'UniformOutput', false);
taxonomyTab = readtable(taxonomyFile, 'readVariableNames', true, 'Delimiter', '\t');
tax_experiment = taxonomyTab(ismember(taxonomyTab.isolate_ID,model_ids),:);
tax_classes = unique(tax_experiment.(2));

% % Export, pathway
% pathways_exported = cell(numel(exported_per_model), 1);
% for i=1:numel(exported_per_model)
%     pathways = map2KEGGPathway(exported_per_model{i}, pathwayFile);
%     pathways_exported{i} = vertcat(pathways{:});
% end
% 
% % shrink the matrix to taxonomic groups
% pathways_exported_group = cell(size(tax_classes));
% for i=1:numel(tax_classes)
%     idx_tax = ismember(tax_experiment.(2), tax_classes{i});
%     pathways_exported_group{i} = vertcat(pathways_exported{idx_tax});
% end
% pathways_exported = pathways_exported_group; clear pathways_exported_group
% 
% % quantify the occurrence of every pathway per model
% all_pathways_exp = unique(vertcat(pathways_exported{:}));
% n_pways_exp = zeros(numel(pathways_exported), numel(all_pathways_exp));
% for i=1:numel(all_pathways_exp)
%     n_pways_exp(:,i) = cellfun(@(x)sum(ismember(x, all_pathways_exp{i})), pathways_exported);
% end

% Export, brite
brite_exported = cell(numel(exported_per_model), 1);
for i=1:numel(exported_per_model)
    brite = map2KEGGBrite(exported_per_model{i}, briteFile);
    brite(cellfun(@isempty, brite)) = {'Other'};
    brite = vertcat(brite{:});
    brite(cellfun(@isempty, brite)) = {'Other'};
    brite_exported{i} = brite;
end

% shrink the matrix to taxonomic groups
brite_exported_group = cell(size(tax_classes));
for i=1:numel(tax_classes)
    idx_tax = ismember(tax_experiment.(2), tax_classes{i});
    brite_exported_group{i} = vertcat(brite_exported{idx_tax});
end
brite_exported = brite_exported_group; clear brite_exported_group

% quantify the occurrence of every brite per model
all_brite_exp = unique(vertcat(brite_exported{:}));
n_brite_exp = zeros(numel(brite_exported), numel(all_brite_exp));
for i=1:numel(all_brite_exp)
    n_brite_exp(:,i) = cellfun(@(x)sum(ismember(x, all_brite_exp{i})), brite_exported);
end

% % Import, pathway
% pathways_imported = cell(numel(imported_per_model), 1);
% for i=1:numel(exported_per_model)
%     pathways = map2KEGGPathway(imported_per_model{i}, pathwayFile);
%     pathways_imported{i} = vertcat(pathways{:});
% end; clear pathways
% 
% % shrink the matrix to taxonomic groups
% pathways_imported_group = cell(size(tax_classes));
% for i=1:numel(tax_classes)
%     idx_tax = ismember(tax_experiment.(2), tax_classes{i});
%     pathways_imported_group{i} = vertcat(pathways_imported{idx_tax});
% end
% pathways_imported = pathways_imported_group; clear pathways_imported_group
% 
% all_pathways_imp = unique(vertcat(pathways_imported{:}));
% n_pways_imp = zeros(numel(pathways_imported), numel(all_pathways_imp));
% for i=1:numel(all_pathways_imp)
%     n_pways_imp(:,i) = cellfun(@(x)sum(ismember(x, all_pathways_exp{i})), pathways_imported);
% end

% Import, brite
brite_imported = cell(numel(imported_per_model), 1);
for i=1:numel(imported_per_model)
    if isempty(imported_per_model{i})
        brite_imported{i} = cell.empty(0,1);
    else
        brite = map2KEGGBrite(imported_per_model{i}, briteFile);
        brite(cellfun(@isempty, brite)) = {{'Other'}};
        brite = vertcat(brite{:});
        brite(cellfun(@isempty, brite)) = {'Other'};
        brite_imported{i} = brite;
    end
end

% shrink the matrix to taxonomic groups
brite_imported_group = cell(size(tax_classes));
for i=1:numel(tax_classes)
    idx_tax = ismember(tax_experiment.(2), tax_classes{i});
    brite_imported_group{i} = vertcat(brite_imported{idx_tax});
end
brite_imported = brite_imported_group; clear brite_imported_group

% quantify the occurrence of every brite per model
all_brite_imp = unique(vertcat(brite_imported{:}));
n_brite_imp = zeros(numel(brite_imported), numel(all_brite_imp));
for i=1:numel(all_brite_imp)
    n_brite_imp(:,i) = cellfun(@(x)sum(ismember(x, all_brite_imp{i})), brite_imported);
end

%% most abundant pathways/brite per unspecific interaction
% pathways_exchange = cell(numel(pathways_exported));
n_groups = numel(tax_classes);
brite_exchange = cell(n_groups);

% get a list of all exported and imported metabolites respectively
all_exp = unique(vertcat(exported_per_model{:}));
all_imp = unique(vertcat(imported_per_model{:}));

% get associated KEGG pathways
% all_exp_pways = map2KEGGPathway(all_exp, pathwayFile);
% all_imp_pways = map2KEGGPathway(all_imp, pathwayFile);

% get associated brite classification
all_exp_brite = map2KEGGBrite(all_exp, briteFile);
all_imp_brite = map2KEGGBrite(all_imp, briteFile);

matrix_exchange_IDs = cell(size(n_brite_imp, 1));

for i=1:numel(brite_exported)
    for j=1:numel(brite_exported)
        % find models belonging to family
        idx_tax_i = ismember(tax_experiment.(2), tax_classes{i});
        idx_tax_j = ismember(tax_experiment.(2), tax_classes{j});
        
        % intersection of export of i and import of j
        exchange =  intersect(...
            vertcat(exported_per_model{idx_tax_i}), vertcat(imported_per_model{idx_tax_j}));
        
%         disp(translateIDs(strtok(exchange, '['), 'met', [], 'MNXref', 'NAMES'))
        if ~isempty(exchange) %&& i~=j
            if exist('all_exp_pways', 'var')
                % associated pathways
                tmp_pways = cellfun(@(x)all_exp_pways(strcmp(all_exp, x)),...
                    exchange);
                tmp_pways = cell(vertcat(tmp_pways{:}));
                
                if ~isempty(tmp_pways)
                    
                    % index associated with the most-abundant pathway
                    [~, tmp_idx] = max(cellfun(@(x)sum(ismember(tmp_pways,x)), tmp_pways));
                    pathways_exchange(i,j) = tmp_pways(tmp_idx);
                else
                    pathways_exchange(i,j) = {''};
                end
            end
            % associated brite
            tmp_brite = cellfun(@(x)all_exp_brite(strcmp(all_exp, x)),...
                exchange);
            tmp_brite = cell(vertcat(tmp_brite{:}));
            tmp_brite(cellfun(@isempty, tmp_brite)) = {'Other'};

            if ~isempty(tmp_brite)
                matrix_exchange_IDs{i,j} = tmp_brite;
                % index associated with the most-abundant pathway/brite
                [~, tmp_idx] = max(cellfun(@(x)sum(ismember(tmp_brite,x)),...
                setdiff(tmp_brite, {'Other'})));
                if ~isempty(tmp_idx)
                    brite_exchange(i,j) = tmp_brite(tmp_idx);
                else
                    brite_exchange(i,j) = {''};
                end
            else
                brite_exchange(i,j) = {''};
            end
        else
            pathways_exchange(i,j) = {''};
            brite_exchange(i,j) = {'Other'};
        end
    end
end
clear exchange tmp_pways tmp_brite tmp_idx


% write to file
% writetable(cell2table(pathways_exchange,...
%     'RowNames', tax_classes, 'VariableNames', tax_classes),...
%     ['/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/graph/',...
%     subFolder, '_pathways_exchanged_',experiment, '.txt'],...
%     'WriteVariableNames', true, 'WriteRowNames', true,...
%     'Delimiter', '\t')

writetable(cell2table(brite_exchange,...
    'RowNames', tax_classes, 'VariableNames', tax_classes),...
    ['/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/graph/',...
    sub_dir, '_brite_exchanged_',experiment, '.txt'],...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter', '\t')

% %% pathway enrichment per model
% alpha = 0.05;
% 
% % scale the occurrences to the number of models per family
% % scale_tax = cellfun(@(x)sum(ismember(tax_experiment.(2), x)), tax_classes);
% % n_pways_exp = bsxfun(@rdivide, n_pways_exp, scale_tax);
% % n_pways_imp = bsxfun(@rdivide, n_pways_imp, scale_tax);
% 
% p_upper_import = zeros(size(n_pways_imp, 1), size(n_pways_imp,2));
% 
% M = sum(sum(n_pways_imp,1));
% 
% for group=1:size(n_pways_imp, 1)
%     for pathway=1:size(n_pways_imp,2)
%         
%         x = n_pways_imp(group, pathway);
%         
%         if x==0
%             p_upper_import(group, pathway) = 1;
%         else
%             K = sum(n_pways_imp(:,pathway));
%             N = sum(n_pways_imp(group, :));
%             
%             p_upper_import(group,pathway) = hygecdf(x,M,K,N, 'upper');
% %             [~, p_upper_import(group,pathway)] = fishertest(...
% %                 [x N-x; K-x M-K-N-x],...
% %                 'Tail', 'right');
%         end
%     end
% end
% 
% % multiple testing correction
% [~, p_adj] = mafdr(reshape(p_upper_import, numel(p_upper_import), 1));
% h_upper_import = p_adj < alpha;
% h_upper_import = reshape(h_upper_import, size(n_pways_imp,1), size(n_pways_imp,2));
% 
% p_upper_export = zeros(size(n_pways_exp, 1), size(n_pways_exp,2));
% 
% M = sum(sum(n_pways_exp,1));
% 
% for group=1:size(n_pways_exp, 1)
%     for pathway=1:size(n_pways_exp,2)
%         
%         x = n_pways_exp(group, pathway);
%         
%         if x==0
%             p_upper_export(group,pathway) = 1;
%         else
%             K = sum(n_pways_exp(:,pathway));
%             N = sum(n_pways_exp(group, :));
%             p_upper_export(group,pathway) = hygecdf(x,M,K,N, 'upper');
% %         [~, p_upper_export(group,pathway)] = fishertest(...
% %             [x N-x; K-x M-K-N+x],...
% %             'Tail', 'right');
%         end
%         
% 
%     end
% end
% 
% % multiple testing correction
% [~, p_adj] = mafdr(reshape(p_upper_export, numel(p_upper_export), 1));
% h_upper_export =  p_adj < alpha;
% h_upper_export = reshape(h_upper_export, size(n_pways_exp,1), size(n_pways_exp,2));
% clear p_fdr
% 
% % find enriched pathways
% enriched_exp = cell(size(h_upper_export,1),1);
% for i=1:size(h_upper_export,1)
%     enriched_exp{i} = all_pathways_exp(h_upper_export(i,:));
% end
% 
% enriched_imp = cell(size(h_upper_import,1),1);
% for i=1:size(h_upper_import,1)
%     enriched_imp{i} = all_pathways_imp(h_upper_import(i,:));
% end

%% pairwise intersection between exported and imported metabolites of each
% two models (i.e. exchanged metabolites)

exchange_mt = zeros(n);
exp_jd_mt = zeros(n);

for i=1:n-1
    for j=i+1:n
        exchange_mt(i,j) = numel(intersect(...
            exported_per_model{i}, imported_per_model{j}));
        exchange_mt(j,i) = numel(intersect(...
            exported_per_model{j}, imported_per_model{i}));
        
        % pairwise Jaccard distance of exported metabolites
        exp_jd_mt(i,j) = numel(intersect(...
            exported_per_model{i}, exported_per_model{j})) / ...
            numel(union(exported_per_model{i}, exported_per_model{j}));
        exp_jd_mt(j,i) = exp_jd_mt(i,j);
    end
end

%% write abundances and exchanged per model to file
% 
% imported_per_group = cell(size(tax_classes));
% exported_per_group = cell(size(tax_classes));
% abundance_per_class = zeros(size(tax_classes));
% for i=1:numel(tax_classes)
%     idx_tax = ismember(tax_experiment.(2), tax_classes{i});
%     imported_per_group{i} = vertcat(imported_per_model{idx_tax});
%     exported_per_group{i} = vertcat(exported_per_model{idx_tax});
%     v(i) = var(cellfun(@numel,exported_per_model(idx_tax)));
%     abundance_per_class(i) = sum(idx_tax);
% end
% 
% writetable(array2table([cellfun(@numel,exported_per_group)./abundance_per_class,...
%     cellfun(@numel,imported_per_group)./abundance_per_class],...
%     'RowNames', tax_classes, 'VariableNames', {'export', 'import'}),...
%     ['/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/',...
%     subFolder, '_exported_vs_abundance_',experiment, '.txt'],...
%     'WriteVariableNames', true, 'WriteRowNames', true,...
%     'Delimiter', '\t')


writetable(array2table([otutab.abundances, cellfun(@numel,exported_per_model)],...
    'RowNames', otutab.Row, 'VariableNames', {'abundance', 'exported'}),...
    ['/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/',...
    sub_dir, '_exported_vs_abundance_',experiment, '.txt'],...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter', '\t')

writetable(array2table([otutab.abundances, cellfun(@numel,imported_per_model)],...
    'RowNames', otutab.Row, 'VariableNames', {'abundance', 'imported'}),...
    ['/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/',...
    sub_dir, '_imported_vs_abundance_',experiment, '.txt'],...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter', '\t')


%% group models together by family
model_ids = cellfun(@(x)strtok(x.id, '_'), GF,...
    'UniformOutput', false);

taxonomyTab = readtable(taxonomyFile, 'readVariableNames', true, 'Delimiter', '\t');
tax_experiment = taxonomyTab(ismember(taxonomyTab.isolate_ID,model_ids),:);
tax_classes = unique(tax_experiment.(2));
n_per_class = cellfun(@(x)sum(contains(tax_experiment.(2), x)), tax_classes);
abundance_per_class = zeros(numel(tax_classes), 1);

for i=1:numel(tax_classes)
    
    % find models that belong to the taxonomic class
    idx_tax = ismember(tax_experiment.(2), tax_classes{i});
    tax_experiment(idx_tax, :) = [];
    % calculate row and column sums for the models belonging to the current
    % class
    row = sum(exchange_mt(idx_tax, :), 1);
    row(idx_tax) = [];
    col = sum(exchange_mt(:, idx_tax), 2);
    col(idx_tax) = [];
    col(end+1) = sum(sum(exchange_mt(idx_tax, idx_tax)));
    
    % remove obsolete rows and columns
    exchange_mt(idx_tax, :) = [];
    exchange_mt(:, idx_tax) = [];
    
    % insert new row and column at the end
    exchange_mt(end+1, :) = row;
    exchange_mt(:, end+1) = col;
    
    % also merge families in exp_jd_mt
    row = mean(exp_jd_mt(idx_tax, :), 1);
    row(idx_tax) = [];
    col = mean(exp_jd_mt(:, idx_tax), 2);
    col(idx_tax) = [];
    col(end+1) = mean(mean(exp_jd_mt(idx_tax, idx_tax)));
    
    % remove obsolete rows and columns
    exp_jd_mt(idx_tax, :) = [];
    exp_jd_mt(:, idx_tax) = [];
    
    % insert new row and column at the end
    exp_jd_mt(end+1, :) = row;
    exp_jd_mt(:, end+1) = col;
    
    % combine the abundances
    abundance_per_class(i) = sum(otutab.abundances(idx_tax));
end


%% scale matrix by abundances of taxonomic classes
exchange_mt = bsxfun(@rdivide, exchange_mt, n_per_class);

abundance_per_class = abundance_per_class ./ n_per_class;

%% write the distance matrix to file
writetable(array2table(exchange_mt, 'VariableNames', tax_classes),...
    ['/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/graph/',...
    [habitat, '_', sub_dir], '_exchanged_metabolites_',experiment, '.txt'],...
    'WriteVariableNames', true,...
    'Delimiter', '\t')



