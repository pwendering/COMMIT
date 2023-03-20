function createMnxrefTranslationTables(mnxref_dir)
%% createMnxrefTranslationTables(mnxref_dir)
% Create tab-delimited files for ID translation using files obtained from
% MetaNetX: https://www.metanetx.org/mnxdoc/mnxref.html
% Input:
%       char mnxref_dir:            directory where the following files are
%                                   located:
%                                   * reac_prop.tsv
%                                   * reac_xref.tsv
%                                   * chem_prop.tsv
%                                   * chem_xref.tsv
% 
% Creates the following files within the given directory:
% * MNXref-met-translation-table.csv
% * MNXref-rxn-translation-table.csv
% * MET_NAMES_FROM_MNXref.csv
% * MNXref_MET_FORMULAE.csv

% check if all required files exist
files = {'chem_prop.tsv', 'chem_xref.tsv', 'reac_prop.tsv', 'reac_xref.tsv'};
for i = 1:numel(files)
    if ~isfile(fullfile(mnxref_dir, files{i}))
        error('File could not be found: %s\n', files{i})
    end
end

%% Metabolites
disp('Generating metabolite translation files')
% read chem_prop.tsv file
% Columns (copied from MetaNetX):
% 1) The identifier of the chemical compound in the MNXref namespace
% 2) The common name of the compound
% 3) The reference compound, i.e. a compound selected in an external 
%    resource that best represent this entry
% 4) The formula of the compound
% 5) The charge of the compound
% 6) The mass of the compound
% 7) The standard InChI of the compound
% 8) The standard InChIKey of the compound
% 9) The SMILES of the compound
chem_prop_table = readtable(fullfile(mnxref_dir, 'chem_prop.tsv'),...
    'FileType', 'text', 'CommentStyle', '#', 'ReadVariableNames', false,...
    'Delimiter', '\t');

% obtain all MNXref metabolite IDs, names, and formulas
disp('Names and Formulas')
mnx_met_ids = chem_prop_table.(1);
mnx_met_names = chem_prop_table.(2);
mnx_met_formulas = chem_prop_table.(4);

if ~ismember('BIOMASS', mnx_met_ids)
    mnx_met_ids = [{'BIOMASS'}; mnx_met_ids];
    mnx_met_names = [{'BIOMASS'}; mnx_met_names];
    mnx_met_formulas = [{''}; mnx_met_formulas];
end

% sort arrays by ID to make sure that the order of there reference IDs
% matches with the order of matched IDs
[mnx_met_ids, id_order] = sort(mnx_met_ids);
mnx_met_names = mnx_met_names(id_order);
mnx_met_formulas = mnx_met_formulas(id_order);

clear chem_prop_table id_order

% read chem_xref.tsv file
% Columns (copied from MetaNetX):
% 1) The identifier of a chemical compound in an external resource
% 2) The corresponding identifier in the MNXref namespace
% 3) The description given by the external resource
chem_xref_table = readtable(fullfile(mnxref_dir, 'chem_xref.tsv'),...
    'FileType', 'text', 'CommentStyle', '#', 'ReadVariableNames', false,...
    'Delimiter', '\t');

% remove rows that are not needed
keep_idx = ...
    startsWith(chem_xref_table.(1), 'kegg.compound:') | ...
    startsWith(chem_xref_table.(1), 'bigg.metabolite:') | ...
    startsWith(chem_xref_table.(1), 'metacyc.compound:') | ...
    startsWith(chem_xref_table.(1), 'seed.compound:') | ...
    startsWith(chem_xref_table.(1), 'CHEBI:');
chem_xref_table = chem_xref_table(keep_idx, :);

% KEGG
disp('MNXref ==> KEGG')
[mnxref, xref] = getXrefForDB(chem_xref_table, 'kegg.compound');
mnx_met_kegg = repmat({''}, size(mnx_met_ids));
mnx_met_kegg(ismember(mnx_met_ids, mnxref)) = xref;
    
% BiGG
disp('MNXref ==> BiGG')
[mnxref, xref] = getXrefForDB(chem_xref_table, 'bigg.metabolite');
mnx_met_bigg = repmat({''}, size(mnx_met_ids));
mnx_met_bigg(ismember(mnx_met_ids, mnxref)) = xref;

% MetaCyc
disp('MNXref ==> MetaCyc')
[mnxref, xref] = getXrefForDB(chem_xref_table, 'metacyc.compound');
mnx_met_metacyc = repmat({''}, size(mnx_met_ids));
mnx_met_metacyc(ismember(mnx_met_ids, mnxref)) = xref;

% ModelSEED
disp('MNXref ==> SEED')
[mnxref, xref] = getXrefForDB(chem_xref_table, 'seed.compound');
mnx_met_modelseed = repmat({''}, size(mnx_met_ids));
mnx_met_modelseed(ismember(mnx_met_ids, mnxref)) = xref;

% ChEBI
disp('MNXref ==> ChEBI')
[mnxref, xref] = getXrefForDB(chem_xref_table, 'CHEBI');
mnx_met_chebi = repmat({''}, size(mnx_met_ids));
mnx_met_chebi(ismember(mnx_met_ids, mnxref)) = xref;

clear chem_xref_table

% combine arrays and write tables
mnx_translation_table = cell2table([mnx_met_ids mnx_met_kegg mnx_met_bigg ...
    mnx_met_metacyc mnx_met_modelseed mnx_met_chebi mnx_met_names],...
    'VariableNames', {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED',...
    'ChEBI', 'NAMES'});
writetable(mnx_translation_table, fullfile(mnxref_dir, 'MNXref-met-translation-table.csv'),...
    'Delimiter', '\t');

mnx_met_names_table = cell2table([mnx_met_ids mnx_met_names]);
writetable(mnx_met_names_table, fullfile(mnxref_dir, 'MET_NAMES_FROM_MNXref.csv'),...
    'Delimiter', '\t', 'WriteVariableNames', false);

mnx_met_formulas_table = cell2table([mnx_met_ids mnx_met_formulas]);
writetable(mnx_met_formulas_table, fullfile(mnxref_dir, 'MNXref_MET_FORMULAE.csv'),...
    'Delimiter', '\t', 'WriteVariableNames', false);
    
clearvars -except mnxref_dir

%% Reactions
disp('Generating reaction translation files')
% read reac_prop.tsv file
% Columns (copied from MetaNetX):
% 1) The identifier of the reaction in the MNXref namespace
% 2) Equation of the reaction in the MNXref namespace (compartmentalized and undirected)
% 3) The original best resource from where this reaction comes from
% 4) The EC(s) associated to this reaction
% 5) Is the equation balanced with respect to elemental composition and charge
% 6) Is this a transport reaction
reac_prop_table = readtable(fullfile(mnxref_dir, 'reac_prop.tsv'),...
    'FileType', 'text', 'CommentStyle', '#', 'ReadVariableNames', false,...
    'Delimiter', '\t');

reac_prop_table = reac_prop_table(~ismember(reac_prop_table.(1), 'EMPTY'),:);

% obtain all MNXref reaction IDs and EC numbers
disp('MNXref ==> EC')
mnx_rxn_ids = reac_prop_table.(1);
mnx_rxn_ec = reac_prop_table.(4);
mnx_rxn_ec = strrep(mnx_rxn_ec, ';', '|');

% sort arrays by ID to make sure that the order of there reference IDs
% matches with the order of matched IDs
[mnx_rxn_ids, id_order] = sort(mnx_rxn_ids);
mnx_rxn_ec = mnx_rxn_ec(id_order);

clear reac_prop_table id_order

% read reac_xref.tsv file
% Columns (copied from MetaNetX):
% 1) The identifier of a biochemical reaction in an external resource
% 2) The corresponding identifier in the MNXref namespace
% 3) The biochemical reaction given by the external resource
reac_xref_table = readtable(fullfile(mnxref_dir, 'reac_xref.tsv'),...
    'FileType', 'text', 'CommentStyle', '#', 'ReadVariableNames', false,...
    'Delimiter', '\t');

% remove rows that are not needed
keep_idx = ...
    startsWith(reac_xref_table.(1), 'kegg.reaction:') | ...
    startsWith(reac_xref_table.(1), 'bigg.reaction:') | ...
    startsWith(reac_xref_table.(1), 'metacyc.reaction:') | ...
    startsWith(reac_xref_table.(1), 'seed.reaction:') | ...
    startsWith(reac_xref_table.(1), 'rhea:');
reac_xref_table = reac_xref_table(keep_idx, :);

reac_xref_table = reac_xref_table(~ismember(reac_xref_table.(2), 'EMPTY'), :);

% KEGG
disp('MNXref ==> KEGG')
[mnxref, xref] = getXrefForDB(reac_xref_table, 'kegg.reaction');
mnx_rxn_kegg = repmat({''}, size(mnx_rxn_ids));
mnx_rxn_kegg(ismember(mnx_rxn_ids, mnxref)) = xref;

% BiGG
disp('MNXref ==> BiGG')
[mnxref, xref] = getXrefForDB(reac_xref_table, 'bigg.reaction');
mnx_rxn_bigg = repmat({''}, size(mnx_rxn_ids));
mnx_rxn_bigg(ismember(mnx_rxn_ids, mnxref)) = xref;

% MetaCyc
disp('MNXref ==> MetaCyc')
[mnxref, xref] = getXrefForDB(reac_xref_table, 'metacyc.reaction');
mnx_rxn_metacyc = repmat({''}, size(mnx_rxn_ids));
mnx_rxn_metacyc(ismember(mnx_rxn_ids, mnxref)) = xref;

% ModelSEED
disp('MNXref ==> SEED')
[mnxref, xref] = getXrefForDB(reac_xref_table, 'seed.reaction');
mnx_rxn_modelseed = repmat({''}, size(mnx_rxn_ids));
mnx_rxn_modelseed(ismember(mnx_rxn_ids, mnxref)) = xref;

% Rhea
disp('MNXref ==> Rhea')
[mnxref, xref] = getXrefForDB(reac_xref_table, 'rhea');
mnx_rxn_rhea = repmat({''}, size(mnx_rxn_ids));
mnx_rxn_rhea(ismember(mnx_rxn_ids, mnxref)) = xref;

clear reac_xref_table

% combine arrays and write table
mnx_translation_table = cell2table([mnx_rxn_ids mnx_rxn_kegg mnx_rxn_bigg ...
    mnx_rxn_metacyc mnx_rxn_modelseed mnx_rxn_rhea mnx_rxn_ec],...
    'VariableNames', {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED',...
    'Rhea', 'EC'});
writetable(mnx_translation_table, fullfile(mnxref_dir, 'MNXref-rxn-translation-table.csv'),...
    'Delimiter', '\t');

%% Test translation
fprintf('Testing translation tables...\n')
passed = testMnxrefTranslationTables(mnxref_dir);
if passed
    fprintf('-- all tests were passed\n')
end

    function [mnxref, xref] = getXrefForDB(xref_table, db_keyword)
    %% [mnxref, xref] = getXrefForDB(xref_table, db_keyword)
    % Obtain a translation from MNXref identifiers to references in another
    % database as specified by the given keyword.
    % Input:
    %   table xref_table:       MNXref chem_xref or reac_xref table
    %   char db_keyword:        keyword to find the lines that correspond
    %                           to the database of interest
    % Output:
    %   cellstr mnxref:         MNXref IDs
    %   cellstr xref:           corresponding IDs in another database
    
    db_idx = startsWith(xref_table.(1), db_keyword);
    mnxref_raw = xref_table.(2)(db_idx);
    xref_raw = erase(xref_table.(1)(db_idx), [db_keyword ':']);
    clear xref_table
    
    % now take the unique, sorted list of MNXref IDs and combine all
    % cross-reference entries
    mnxref = unique(mnxref_raw);
    xref = cellfun(@(x)cellstr(strjoin(xref_raw(ismember(mnxref_raw, x)), '|')),...
        mnxref);
    
    end

end