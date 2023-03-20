function passed = testMnxrefTranslationTables(mnxref_dir)
%% passed = testMnxrefTranslationTables(mnxref_dir)
% Input:
%  char mnxref_dir:            directory where the following files are
%                              located:
%                              * MNXref-met-translation-table.csv
%                              * MNXref-rxn-translation-table.csv
%                              * MET_NAMES_FROM_MNXref.csv
%                              * MNXref_MET_FORMULAE.csv
% Ouput:
%  logical passed:             true if ID matching test was passing

passed = false;

%% Metabolites
met_format_string = '%s\t%s\t%s\t%s\t%s\t%s\t%s';

mnxref_test_ids     = {'MNXM1', 'MNXM10', 'MNXM3', 'MNXM1137670', 'WATER',...
    'MNXM1062', 'MNXM1106852', 'MNXM740811', 'MNXM9132', 'MNXM1107753',...
    'MNXM124', 'MNXM1102167'};
kegg_test_ids       = {'C00080', 'C00004', 'C00002', 'C00031', 'C00001',...
    'C01250', 'C16535', 'C11378', 'C05463', 'C00158', 'C00315', 'C00021'};
bigg_test_ids       = {'h', 'nadh', 'atp', 'glc__D', 'h2o', 'acg5sa',...
    'M02613', 'q10', 'tdechola', 'cit', 'spmd', 'ahcys'};
metacyc_test_ids    = {'PROTON', 'NADH', 'ATP', 'Glucopyranose', 'WATER',...
    'CPD-469', 'CPD-18572', 'UBIQUINONE-10', 'CPD-12457', 'CIT',...
    'SPERMIDINE', 'ADENOSYL-HOMO-CYS'};
modelseed_test_ids  = {'cpd00067', 'cpd00004', 'cpd00002', 'cpd00027',...
    'cpd00001', 'cpd00918', 'cpd16349', 'cpd08232', 'cpd03244', 'cpd00137',...
    'cpd00264', 'cpd00019'};
chebi_test_ids      = {'15378', '57945', '30616', '4167', '33813', '29123',...
    '78796', '46245', '36261', '16947', '57834', '57856'};
name_test_ids       = {'H(+)', 'NADH', 'ATP', 'D-glucose', 'H2O',...
    'N-acetyl-L-glutamate 5-semialdehyde', 'nonadecanoate', 'ubiquinone-10',...
    'taurodeoxycholate', 'citrate', 'spermidine', 'S-adenosyl-L-homocysteine'};

test_db_idx = {kegg_test_ids, bigg_test_ids, metacyc_test_ids,...
    modelseed_test_ids, chebi_test_ids, name_test_ids};

mnx_met_t_file = fullfile(mnxref_dir, 'MNXref-met-translation-table.csv');

% Specify range and delimiter
mnxref_met_tab = readtable(mnx_met_t_file,...
    'Format', met_format_string);

% check, if the first lines were parsed properly
if ~all(ismember({'BIOMASS', 'MNXM01', 'MNXM02'}, mnxref_met_tab.MNXref))
    warning(['The metabolite translation table was not read properly, '...
        'IDs BIOMASS, MNXM01, and MNXM02 are missing'])
    return
end

for i = 1:numel(mnxref_test_ids)
    idx = ismember(mnxref_met_tab.MNXref, mnxref_test_ids(i));
    if sum(idx) == 1
        for j = 1:numel(test_db_idx)
            % sometimes one MetaNetX ID matches to multiple database IDs,
            % only check if the specified test ID is one of them
            db_id = test_db_idx{j}{i};
            id_match = ismember(db_id,...
                strsplit(mnxref_met_tab.(j+1){idx}, '|'));
            if ~id_match
                warning('Expected cross-reference %s was not found for ID %s',...
                    db_id, mnxref_test_ids{i})
                return
            end
        end
    else
        warning('ID %s was not found or has multiple occurences',...
            mnxref_test_ids{1})
        return
    end
end

clearvars -except mnxref_dir passed

%% Reactions
rxn_format_string = '%s\t%s\t%s\t%s\t%s\t%s\t%s';

mnxref_test_ids     = {'MNXR100024', 'MNXR100450', 'MNXR100492', 'MNXR100947',...
    'MNXR104272', 'MNXR109319', 'MNXR112092', 'MNXR147173', 'MNXR189729',...
    'MNXR97932', 'MNXR99896'};
kegg_test_ids       = {'R00253', 'R00497', 'R07172', 'R01960', 'R00610',...
    'R05102', 'R08558', 'R01176', 'R01186', 'R00658', 'R03313'};
bigg_test_ids       = {'GALh', 'GTHS', 'H2O2syn', 'KYN3OX', 'SARCOXp',...
    'CYTP450Rh', 'CHOLD3', 'HMR_0156', 'MI4PP', 'ENO_h', 'G5DHx'};
metacyc_test_ids    = {'GLUTAMINESYN-RXN', 'GLUTATHIONE-SYN-RXN', '',...
    'KYNURENINE-3-MONOOXYGENASE-RXN', 'SARCOX-RXN', '', 'RXN-12582',...
    'RXN-19193', 'RXN-10952', '2PGADEHYDRAT-RXN', 'GLUTSEMIALDEHYDROG-RXN'};
modelseed_test_ids  = {'rxn33952', 'rxn00351', 'rxn04958', 'rxn13171',...
    'rxn31075', 'rxn06853', 'rxn12191', 'rxn00873', 'rxn00882', 'rxn00459',...
    'rxn02373'};
rhea_test_ids      = {'16169', '13557', '11260', '20545', '13313', '', '',...
    '46172', '30735', '10164', '19541'};
ec_test_ids       = {'6.3.1.2', '6.3.2.3', '1.6.3.1', '1.14.13.9', '1.5.3.1',...
    '1.6.2.4', '1.1.1.1', '6.2.1.2', '3.1.3.25', '4.2.1.11', '1.2.1.41'};

test_db_idx = {kegg_test_ids, bigg_test_ids, metacyc_test_ids,...
    modelseed_test_ids, rhea_test_ids, ec_test_ids};

mnx_rxn_t_file = fullfile(mnxref_dir, 'MNXref-rxn-translation-table.csv');

% Specify range and delimiter
mnxref_rxn_tab = readtable(mnx_rxn_t_file,...
    'Format', rxn_format_string);

% check, if the first lines were parsed properly
if ~all(ismember({'MNXR01', 'MNXR02', 'MNXR03'}, mnxref_rxn_tab.MNXref))
    warning(['The reaction translation table was not read properly, '...
        'IDs MNXR01, MNXR02, and MNXR03 are missing'])
    return
end

if ismember('EMPTY', mnxref_rxn_tab.MNXref)
    warning('MetaNetX reaction ID ''EMPTY'' detect, which is not allowed.')
    return
end

for i = 1:numel(mnxref_test_ids)
    idx = ismember(mnxref_rxn_tab.MNXref, mnxref_test_ids(i));
    if sum(idx) == 1
        for j = 1:numel(test_db_idx)
            % sometimes one MetaNetX ID matches to multiple database IDs,
            % only check if the specified test ID is one of them
            db_id = test_db_idx{j}{i};
            id_match = ismember(db_id,...
                strsplit(mnxref_rxn_tab.(j+1){idx}, '|'));
            if ~id_match
                warning('Expected cross-reference %s was not found for ID %s',...
                    db_id, mnxref_test_ids{i})
                return
            end
        end
    else
        warning('ID %s was not found or has multiple occurences',...
            mnxref_test_ids{1})
        return
    end
end

passed = true;

end

