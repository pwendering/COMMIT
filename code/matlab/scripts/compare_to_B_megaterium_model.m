%% compare to Bacillus megaterium metabolic model as reference
warning 'off' 
wd = '/stud/wendering/Masterthesis/DATA/reference-models';

tablesDir = '/stud/wendering/Masterthesis/DATA/tables';

plotDir = '/stud/wendering/Masterthesis/FIGURES/Reference-models';
modelDir = fullfile(wd, 'models_Zou_2013');
phyloFile = '/stud/wendering/Masterthesis/DATA/genomes/B_megaterium/16S-seqs.nw.dist.txt';
ecTransFile = fullfile(tablesDir, 'corrected-EC-numbers.csv');
ecTranslationTable = readtable(ecTransFile, 'ReadVariableNames', false);
%% Pre-processing

% Reference: Zou et al., 2013, J. Biotechnol.
% uploaded to MetaNetX
% reference = readCbModel(fullfile(wd, 'models_EBI', 'Bacillus-megaterium-WSH002-BMID000000142083_translated.xml'));
% reference.rxns = reference.rxnMetaNetXID;
% reference.mets = reference.metMetaNetXID;

reference = readCbModel(fullfile(modelDir, 'Bacillus-megaterium-WSH-200-iMZ1055.xml'));
% reference = removeRxns(reference, reference.rxns(~(cellfun('isempty', regexp(reference.rxns, '_(in)|(out)|(bm)$')))));
reference.description = 'B-megaterium-ref';

% Genes
% reference.genes = regexprep(reference.genes, '_i$', '');

% Reactions
rxns_ref = {};
% rxns_ref = regexp(reference.rxnKEGGID, '; ', 'split');
% rxns_ref = [rxns_ref{:}];
% rxns_ref = translateIDs(rxns_ref, 'rxn', translationDB.rxnTab,...
%     'KEGG', 'MNXref');

% Metabolites
mets_ref = {};
% mets_ref = regexp(reference.metNotes, 'KEGG:(\ [CDG]\d{5},*)+', 'match');
% mets_ref = [mets_ref{:}];
% mets_ref = regexp(mets_ref, '[CG]\d{5}', 'match');
% mets_ref = [mets_ref{:}];
% mets_ref = translateIDs(mets_ref, 'met', translationDB.metTab,...
%     'KEGG', 'MNXref');


% EC numbers
ec_ref = regexp(reference.rxnECNumbers, '/', 'split');
ec_ref = [ec_ref{:}];
ec_ref = correctEC(ec_ref, ecTranslationTable);
ec_ref = regexp(ec_ref, '|', 'split');
ec_ref = [ec_ref{:}];
ec_ref = ec_ref(~cellfun('isempty', ec_ref));
% translate EC number to MNXref reactions
rxns_ref = translateIDs(ec_ref, 'rxn', [], 'EC', 'MNXref', false);

adj_rxns = 1 - sum(findExcRxns(reference)) / numel(reference.rxns);
% corresponding draft models:
% same species: Leaf75, Soil531
% same genus:
%   Root920, Root11, Root131, Root239
%   Soil768D1, Soil745
%   Leaf13, Leaf49, Leaf406

%% Initialization of count variables
n_models = 432;
TP = zeros(n_models, 1); % intersection of model and reference
FP = zeros(n_models, 1); % contained in the model but not in the reference
FN = zeros(n_models, 1); % contained in the reference but not in the model

model_ids = {};
% same species
species_ids = {'Leaf75', 'Soil531'};
% same genus
genus_ids = {'Leaf13', 'Leaf49', 'Leaf406',...
    'Root920', 'Root11', 'Root131', 'Root239',...
    'Soil768D1', 'Soil745'};

p_rxn = zeros(n_models, 1);
% p_met = zeros(n_models, 1);
p_ec = zeros(n_models, 1);
svd_dist = zeros(n_models, 1);
p_species_genes = zeros(size(species_ids));
p_genus_genes = zeros(size(genus_ids));

% parpool(2);

c = 0;
for habitat = {'Leaf', 'Root', 'Soil'}
    
    % consensus
     spec = '';
     load(['/stud/wendering/Masterthesis/DATA/Consensus_models/',...
         char(habitat), '_consensus_models.mat'])
    % KBase draft
    %spec = 'KBase_';
    %load(fullfile('/stud/wendering/Masterthesis/DATA/models_KBase',...
    %    char(habitat), [char(habitat), '_models_genes_translated.mat']))
    % RAVEN 2.0 draft
    %spec = 'RAVEN_';
    %load(fullfile('/stud/wendering/Masterthesis/DATA/models_RAVEN/HMMer10E-50/',...
    %   [char(habitat), '_models_no_medium_no_biomass']))
    % CarveMe draft
    %spec = 'CarveMe_';
    %load(fullfile('/stud/wendering/Masterthesis/DATA/models_CarveMe',...
    %    char(habitat), [char(habitat), '_models_no_medium_no_biomass']))
    % AuReMe draft
    %spec = 'AuReMe_';
    %load(fullfile('/stud/wendering/Masterthesis/DATA/models_AuReMe',...
    %    char(habitat), [char(habitat), '_models_genes_translated']))
    if ~exist('merged_models', 'var')
	    merged_models = models;
    end
    merged_models = reshape(merged_models, numel(merged_models), 1);
    current_ids = cellfun(@(x)strtok(x.id, '_'), merged_models, 'un', 0);
    model_ids = [model_ids; current_ids];
    disp(habitat)

    for i=1:numel(merged_models)
        
        c = c + 1;
        
        rxns_draft = strtok(merged_models{i}.rxns, '_');
        mets_draft = strtok(merged_models{i}.mets, '[');
        ec_draft = regexp(merged_models{i}.EC, '\|', 'split');
        ec_draft = [ec_draft{:}];
        
        I_rxn = intersect(rxns_draft, rxns_ref);
        I_met = intersect(mets_draft, mets_ref);
        I_ec = intersect(ec_draft, ec_ref);
        
        TP(c) = numel(I_rxn) + ...
            numel(I_met) + ...
            numel(I_ec);
        FN(c) = numel(setdiff(rxns_ref, rxns_draft)) + ...
            numel(setdiff(mets_ref, mets_draft)) + ...
            numel(setdiff(ec_ref, ec_draft));
        FP(c) = numel(setdiff(rxns_draft, rxns_ref)) + ...
            numel(setdiff(mets_draft, mets_ref)) + ...
            numel(setdiff(ec_draft, ec_ref));
        
        p_rxn(c) = numel(I_rxn) / ( numel(rxns_ref) * adj_rxns );
%         p_met(c) = numel(I_met) / numel(mets_ref);
        p_ec(c) = numel(I_ec) / ( numel(ec_ref) * adj_rxns );
        
    end
    
    % SVD distance
    %svd_mt = svdDistance([merged_models', reference]);
    %svd_dist(c-numel(merged_models)+1:c) = svd_mt(1:end-1, end);
        
    % Map gene IDs:
    
    % same species
    for j = find(contains(species_ids, habitat))
        filename = strcat('/stud/wendering/Masterthesis/DATA/genomes/B_megaterium/B-megaterium-',...
            species_ids{j}, '.mapping');
        model = translateGeneIDs(reference, filename);
        
        p_species_genes(j) = numel(intersect(model.genes,...
            merged_models{ismember(current_ids, species_ids{j})}.genes)) / ...
            numel(model.genes); clear model
    end
    
    % same genus
    for j = find(contains(genus_ids, habitat))
        filename = strcat('/stud/wendering/Masterthesis/DATA/genomes/B_megaterium/B-megaterium-',...
            genus_ids{j}, '.mapping');
        model = translateGeneIDs(reference, filename);
        
        p_genus_genes(j) = numel(intersect(model.genes,...
            merged_models{ismember(current_ids, genus_ids{j})}.genes)) / ...
            numel(model.genes); clear model
    end
clear merged_models 
end

clear current_ids svd_mt

idx_species = ismember(model_ids, species_ids);
idx_genus = ismember(model_ids, genus_ids);

p_genus_rxn = p_rxn(idx_genus);
% p_genus_met = p_met(idx_genus);
p_genus_ec = p_ec(idx_genus);

p_species_rxn = p_rxn(idx_species);
% p_species_met = p_met(idx_species);
p_species_ec = p_ec(idx_species);

% Quality measures
prec = TP ./ (TP + FP);
sens = TP ./ (TP + FN);
F1 = 2*TP ./ (2*TP + FN + FP);

prec_species = prec(idx_species);
sens_species = sens(idx_species);
F1_species = F1(idx_species);

prec_genus = prec(idx_genus);
sens_genus = sens(idx_genus);
F1_genus = F1(idx_genus);


%svd_species = svd_dist(idx_species);
%svd_genus = svd_dist(idx_genus);

% Phylogenetic distance:
phylo_dist = readtable(phyloFile,...
    'ReadVariableNames', true,...
    'ReadRowNames', true);
row_names = phylo_dist.Properties.RowNames;
phylo_dist = phylo_dist.BMWSH_16S01;
phylo_dist(contains(row_names, 'BMWSH_16S01')) = [];
row_names(contains(row_names, 'BMWSH_16S01')) = [];
re_order = cellfun(@(x)find(strcmp(x, row_names)), model_ids);
relation = 1 - phylo_dist(re_order); clear re_order row_names phylo_dist

relation_species = relation(idx_species);
relation_genus = relation(idx_genus);

%% Plots

% Boxplot
% boxplot([p_rxn, p_ec, prec, sens, F1],...
%     {'reactions', 'E.C. numbers', 'precision', 'sensitivity', 'F1-score'},...
%     'Symbol', 'ko', 'Colors', 'k')
% title('Distribution of reference model features contained in the consensus models')
% hold on
% % scatter(ones(size(p_rxn)).*(1+(rand(size(p_rxn))-0.5)/5), p_rxn, 'r', 'filled')
% scatter(ones(size(p_genus_rxn)), p_genus_rxn, 'b', 'filled')
% scatter(ones(size(p_species_rxn)), p_species_rxn, 'r', 'filled')
% 
% scatter(repmat(2, size(p_genus_met)), p_genus_met, 'b', 'filled')
% scatter(repmat(2, size(p_species_met)), p_species_met, 'r', 'filled')
% 
% scatter(repmat(3, size(p_genus_ec)), p_genus_ec, 'b', 'filled')
% scatter(repmat(3, size(p_species_ec)), p_species_ec, 'r', 'filled')
% 
% scatter(repmat(4, size(prec_genus)), prec_genus, 'b', 'filled')
% scatter(repmat(4, size(prec_species)), prec_species, 'r', 'filled')
% 
% scatter(repmat(5, size(sens_genus)), sens_genus, 'b', 'filled')
% scatter(repmat(5, size(sens_species)), sens_species, 'r', 'filled')
% 
% scatter(repmat(6, size(F1_genus)), F1_genus, 'b', 'filled')
% scatter(repmat(6, size(F1_species)), F1_species, 'r', 'filled')




% identity vs. relatedness
%subplot(2,2,1)
%scatter(relation, p_ec)
%hold on
%scatter(relation_genus, p_genus_ec, 'b', 'filled')
%scatter(relation_species, p_species_ec, 'r', 'filled')
%title('Proportion of E.C. numbers of reference model contained in draft vs. sequence similarity')
%xlabel('1 - phylogenetic distance')
%ylabel('Proportion of contained E.C. numbers')

% subplot(2,2,1)
% scatter(relation, p_rxn)
% hold on
% scatter(relation_species, p_species_rxn, 'r', 'filled')
% scatter(relation_genus, p_genus_rxn, 'b', 'filled')
% title('Proportion of reactions of reference model contained in draft vs. sequence similarity')
% xlabel('1 - phylogenetic distance')
% ylabel('Proportion of contained reactions')

% subplot(2,2,1)
% scatter(relation, svd_dist)
% hold on
% scatter(relation_genus, svd_genus, 'b', 'filled')
% scatter(relation_species, svd_species, 'r', 'filled')
% title('SVD distance  vs. sequence similarity')
% xlabel('1 - phylogenetic distance')
% ylabel('SVD distance')

%subplot(2,2,2)
%scatter(relation, sens)
%hold on
%scatter(relation_genus, sens_genus, 'b', 'filled')
%scatter(relation_species, sens_species, 'r', 'filled')
%title('Sensitivity vs. sequence similarity')
%xlabel('1 - phylogenetic distance')
%ylabel('sensitivity')

%subplot(2,2,3)
%scatter(relation, prec)
%hold on
%scatter(relation_genus, prec_genus, 'b', 'filled')
%scatter(relation_species, prec_species, 'r', 'filled')
%title('Precision vs. sequence similarity')
%xlabel('1 - phylogenetic distance')
%ylabel('precision')

%subplot(2,2,4)
%scatter(relation, F1)
%hold on
%scatter(relation_genus, F1_genus, 'b', 'filled')
%scatter(relation_species, F1_species, 'r', 'filled')
%title('F1-score vs. sequence similarity')
%xlabel('1 - phylogenetic distance')
%ylabel('F1-score')

%hold off
%scatter(relation_genus, p_genus_genes)
%hold on
%scatter(relation_species, p_species_genes, 'r', 'filled')
%title('Percentage reference model genes the are contained in the merged models')
%xlabel('1 - phylogenetic distance')
%ylabel('Percenage of genes contained')

%saveas(gcf, fullfile(plotDir, 'B_megaterium_genes_vs_phylo.png'))

% saveas(gcf, fullfile(plotDir, 'B_megaterium_dist_vs_phylo.png'))
% save(fullfile(modelDir, 'model_comparison_B_megaterium'),...
%     'sens', 'prec', 'p_ec', 'p_rxn', 'F1', 'svd_dist', 'relation', 'TP',...
%     'FP', 'FN', 'model_ids')



%load(fullfile(modelDir, 'model_comparison_B_megaterium'))
writetable(array2table([p_ec, sens, prec, relation, idx_genus, idx_species],...
    'VariableNames', {'p_EC', 'sensitivity', 'precision', 'relation', 'genus', 'species'},...
    'RowNames', model_ids),...
    ['/stud/wendering/Masterthesis/FIGURES/Comparison-to-reference-models/', spec, 'B_megaterium.txt'],...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter', '\t')






