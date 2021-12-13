%% compare to Bacillus megaterium metabolic model as reference
options

figOutDir = 'figures/Comparison-to-reference-models';
modelDir = fullfile('data', 'reference-models', 'models_Zou_2013');
phyloFile = 'data/genomes/B_megaterium/16S-seqs.nw.dist.txt';
ecTranslationTable = readtable('data/tables/corrected-EC-numbers.csv',...
    'ReadVariableNames', false);

%% Pre-processing

% Reference: Zou et al., 2013, J. Biotechnol.
% read reference model SBML file
reference = readCbModel(fullfile(modelDir, 'Bacillus-megaterium-WSH-200-iMZ1055.xml'));
reference.description = 'B-megaterium-ref';

% Metabolites
mets_ref = table2cell(readtable(fullfile(modelDir,'kegg_met_ids.csv')));
mets_ref = translateIDs(mets_ref, 'met', [], 'KEGG', 'MNXref', false);

% proportion of metabolites that can be compared
adj_met = numel(unique(mets_ref)) / numel(unique(strtok(reference.mets, '[')));

% EC numbers
ec_ref = regexp(reference.rxnECNumbers, '/', 'split');
ec_ref = [ec_ref{:}];
ec_ref = correctEC(ec_ref, ecTranslationTable);
ec_ref = regexp(ec_ref, '|', 'split');
ec_ref = [ec_ref{:}];
ec_ref = ec_ref(~cellfun('isempty', ec_ref));

adj_ec = numel(ec_ref) / numel(reference.rxnECNumbers);

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
p_ec = zeros(n_models, 1);
svd_dist = zeros(n_models, 1);
p_species_genes = zeros(size(species_ids));
p_genus_genes = zeros(size(genus_ids));

c = 0;
for habitat = {'Leaf', 'Root', 'Soil'}
    
    % consensus
    spec = '';
    load(fullfile('data/models/consensus', [char(habitat),...
       '_consensus_models.mat']))
    % KBase draft
    % spec = 'KBase_';
    % load(fullfile('data/models/kbase', char(habitat),...
    %    [char(habitat), '_models_genes_translated.mat']))
    % RAVEN 2.0 draft
    % spec = 'RAVEN_';
    % load(fullfile('data/models/raven/HMMer10E-50/',...
    %    [char(habitat), '_models_no_medium_no_biomass']))
    % CarveMe draft
    % spec = 'CarveMe_';
    % load(fullfile('data/models/carveme',...
    %     char(habitat), [char(habitat), '_models_no_medium_no_biomass']))
    % AuReMe draft
    % spec = 'AuReMe_';
    % load(fullfile('data/models/aureme',...
    %     char(habitat), [char(habitat), '_models_genes_translated']))
    
    if ~exist('merged_models', 'var')
        merged_models = models;
    end
    
    merged_models = reshape(merged_models, numel(merged_models), 1);
    current_ids = cellfun(@(x)strtok(x.id, '_'), merged_models, 'un', 0);
    model_ids = [model_ids; current_ids];
    
    disp(habitat)
    
    for i=1:numel(merged_models)
        
        c = c + 1;
        
        mets_draft = strtok(merged_models{i}.mets, '[');
        ec_draft = regexp(merged_models{i}.EC, '\|', 'split');
        ec_draft = [ec_draft{:}];
        
        I_met = intersect(mets_draft, mets_ref);
        I_ec = intersect(ec_draft, ec_ref);
        
        % true positives
        TP(c) = numel(I_ec) + numel(I_met);
        % false negatives
        FN(c) = numel(setdiff(ec_ref, ec_draft)) + numel(setdiff(mets_ref, mets_draft));
        % false positives
        FP(c) = numel(setdiff(ec_draft, ec_ref)) + numel(setdiff(mets_draft, mets_ref));
        
        p_ec(c) = numel(I_ec) / ( numel(ec_ref) * adj_ec );
    end
    %{
    % Map gene IDs:
    % same species
    for j = find(contains(species_ids, habitat))
        filename = strcat('/data/genomes/B_megaterium/B-megaterium-',...
            species_ids{j}, '.mapping');
        model = translateGeneIDs(reference, filename);
        
        p_species_genes(j) = numel(intersect(model.genes,...
            merged_models{ismember(current_ids, species_ids{j})}.genes)) / ...
            numel(model.genes); clear model
    end
    
    % same genus
    for j = find(contains(genus_ids, habitat))
        filename = strcat('/data/genomes/B_megaterium/B-megaterium-',...
            genus_ids{j}, '.mapping');
        model = translateGeneIDs(reference, filename);
        
        p_genus_genes(j) = numel(intersect(model.genes,...
            merged_models{ismember(current_ids, genus_ids{j})}.genes)) / ...
            numel(model.genes); clear model
    end
    %}
    clear merged_models
end

clear current_ids svd_mt

idx_species = ismember(model_ids, species_ids);
idx_genus = ismember(model_ids, genus_ids);

% Quality measures
precision = TP ./ (TP + FP);
sensitivity = TP ./ (TP + FN);
F1 = 2*TP ./ (2*TP + FN + FP);

prec_species = precision(idx_species);
sens_species = sensitivity(idx_species);
F1_species = F1(idx_species);

prec_genus = precision(idx_genus);
sens_genus = sensitivity(idx_genus);
F1_genus = F1(idx_genus);

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

% write results to file
writetable(array2table([p_ec, sensitivity, precision, relation, idx_genus, idx_species],...
    'VariableNames', {'p_EC', 'sensitivity', 'precision', 'relation', 'genus', 'species'},...
    'RowNames', model_ids),...
    fullfile(figOutDir, [spec, 'B_megaterium.txt']),...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter', '\t')

