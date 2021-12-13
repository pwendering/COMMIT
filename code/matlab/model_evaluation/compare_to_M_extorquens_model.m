% compare to Methylobacterium extorquens metabolic model as reference
options


figOutDir = 'figures/Comparison-to-reference-models';
modelDir = 'data/reference-models/models_Peyraud_BMC_2011';
phyloFile = 'data/genomes/M_extorquens/16S-seqs.nw.dist.txt';
ecTranslationTable = readtable('data/tables/corrected-EC-numbers.csv',...
    'ReadVariableNames', false);

%% Methylobacterium extorquens
% Reference: Peyraud et al., 2011, BMC Syst. Biol.
% read reference model from SBML file
reference = readCbModel(fullfile(modelDir, 'Methylobacterium-extorquens-iRP911.xml'));
reference = removeRxns(reference, reference.rxns(~(cellfun('isempty', regexp(reference.rxns, 'EX_')))));
reference.description = 'M_extorquens-ref';

% Combination of KEGG and MetaCyc IDs from excel sheet (suppl.)
mets_ref = table2cell(readtable(fullfile(modelDir, 'iRP911_mets_KEGG_MCyc.csv'),...
    'ReadVariableNames', false));
mets_ref = translateIDs(mets_ref, 'met', [], 'KEGG', 'MNXref');
mets_ref = translateIDs(mets_ref, 'met', [], 'MetaCyc', 'MNXref');

% proportion of metabolites that can be compared
adj_met = numel(mets_ref) / numel(unique(strtok(reference.mets, '[')));

ec_ref = table2cell(readtable(fullfile(modelDir, 'iRP911_ec_uniq.csv'),...
    'ReadVariableNames', false));
ec_ref = regexp(ec_ref, '\d+\.\d+\.\d+\.\d+', 'match');
ec_ref = correctEC(ec_ref, ecTranslationTable);
ec_ref = [ec_ref{:}];
ec_ref = regexp(ec_ref, '|', 'split');
ec_ref = [ec_ref{:}];
ec_ref = ec_ref(~cellfun('isempty', ec_ref));

adj_ec = numel(ec_ref) / numel(reference.rxns);

reference.genes = table2cell(readtable(fullfile(modelDir, 'iRP911_genes_uniq.csv'),...
    'ReadVariableNames', false, 'Delimiter', '\t'));


% same species:
% Leaf90, Leaf92, Leaf119, Leaf121, Leaf122

% same genus:
% Leaf85, Leaf86, Leaf87, Leaf88, Leaf89, Leaf91, Leaf93, Leaf94, Leaf99,
% Leaf100, Leaf102, Leaf104, Leaf106, Leaf108, Leaf111, Leaf113, Leaf117,
% Leaf123, Leaf125, Leaf344, Leaf361, Leaf399, Leaf456, Leaf465, Leaf466,
% Leaf469, Root483D1

%% Initialization of count variables
n_models = 432;
TP = zeros(n_models, 1); % intersection of model and reference
FP = zeros(n_models, 1); % contained in the model but not in the reference
FN = zeros(n_models, 1); % contained in the reference but not in the model
c = 0;
model_ids = {};

% same species
species_ids = {'Leaf90', 'Leaf92', 'Leaf119', 'Leaf121', 'Leaf122'};

% same genus
genus_ids = {'Leaf85', 'Leaf86', 'Leaf87', 'Leaf88', 'Leaf89', 'Leaf91',...
    'Leaf93', 'Leaf94', 'Leaf99', 'Leaf100', 'Leaf102', 'Leaf104', 'Leaf106',...
    'Leaf108', 'Leaf111', 'Leaf113', 'Leaf117', 'Leaf123', 'Leaf125', 'Leaf344',...
    'Leaf361', 'Leaf399', 'Leaf456', 'Leaf465', 'Leaf466', 'Leaf469', 'Root483D1'};

p_rxn = zeros(n_models, 1);
p_met = zeros(n_models, 1);
p_ec = zeros(n_models, 1);
p_species_genes = zeros(size(species_ids));
p_genus_genes = zeros(size(genus_ids));

for habitat = {'Leaf', 'Root', 'Soil'}
    
    % consensus
    spec = '';
    load(fullfile('data/models/consensus', [char(habitat),...
        '_consensus_models.mat']))
    % KBase draft
    % spec = 'KBase_';
    % load(fullfile('data/models/kbase', char(habitat),...
    %  [char(habitat), '_models_genes_translated.mat']))
    % RAVEN 2.0 draft
    % spec = 'RAVEN_';
    % load(fullfile('data/models/raven/HMMer10E-50/',...
    %  [char(habitat), '_models_no_medium_no_biomass']))
    % CarveMe draft
    % spec = 'CarveMe_';
    % load(fullfile('data/models/carveme',...
    %   char(habitat), [char(habitat), '_models_no_medium_no_biomass']))
    % AuReMe draft
    % spec = 'AuReMe_';
    % load(fullfile('data/models/aureme',...
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
        
        mets_draft = strtok(merged_models{i}.mets, '[');
        ec_draft = regexp(merged_models{i}.EC, '\|', 'split');
        ec_draft = [ec_draft{:}];
        
        I_met = intersect(mets_draft, mets_ref);
        I_ec = intersect(ec_draft, ec_ref);
        
        % true positives
        TP(c) = numel(I_met) + numel(I_ec);
        % false negatives
        FN(c) = numel(setdiff(mets_ref, mets_draft)) + numel(setdiff(ec_ref, ec_draft));
        % false positives
        FP(c) = numel(setdiff(mets_draft, mets_ref)) + numel(setdiff(ec_draft, ec_ref));
        
        p_ec(c) = numel(I_ec) / (numel(ec_ref)*adj_ec);
        
    end
    
    %{
     % Map gene IDs:
     
     % same species
     for j = find(contains(species_ids, habitat))
         filename = strcat('/data/genomes/M_extorquens/M-extorquens-AM1-',...
             species_ids{j}, '.mapping');
         model = translateGeneIDs(reference, filename);
         
         p_species_genes(j) = numel(intersect(model.genes,...
             merged_models{ismember(current_ids, species_ids{j})}.genes)) / ...
             numel(model.genes); clear model
     end
     
     % same genus
     for j = find(contains(genus_ids, habitat))
         filename = strcat('data/genomes/M_extorquens/M-extorquens-AM1-',...
             genus_ids{j}, '.mapping');
         model = translateGeneIDs(reference, filename);
         
         p_genus_genes(j) = numel(intersect(model.genes,...
             merged_models{ismember(current_ids, genus_ids{j})}.genes)) / ...
             numel(model.genes); clear model
     end
    %}
    clear merged_models
end

clear current_ids

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
phylo_dist = phylo_dist.META1p16S5;
phylo_dist(contains(row_names, 'META1p16S5')) = [];
row_names(contains(row_names, 'META1p16S5')) = [];
re_order = cellfun(@(x)find(strcmp(x, row_names)), model_ids);
relation = 1 - phylo_dist(re_order); clear re_order row_names phylo_dist

relation_species = relation(idx_species);
relation_genus = relation(idx_genus);

%% write results to file
writetable(array2table([p_ec, sensitivity, precision, relation, idx_genus, idx_species],...
    'VariableNames', {'p_EC', 'sensitivity', 'precision', 'relation', 'genus', 'species'},...
    'RowNames', model_ids),...
    [figOutDir, filesep, spec, 'M_extorquens.txt'],...
    'WriteVariableNames', true, 'WriteRowNames', true,...
    'Delimiter', '\t')


