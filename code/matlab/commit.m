% load options
options

%% RAVEN draft model reconstruction
model_generation_RAVEN

% KBase, CarveMe, and AuReMe/Pathway Tools draft models were created outside of
% MATLAB and are available upon request

%% Convert draft models to a common format that can be used for evaluation

% % Convert and translate the draft models to MNXref namespace
modelConversion_RAVEN
modelConversion_KBase
modelConversion_CarveMe
modelConversion_AuReMe

%% Evaluate the draft models from all approaches

% Highly-abundant metabolites (connect reactions/EC numbers
% on the basis of connecting metabolites in S)
% E. coli: Highly connected enzymes tend to be highly abundant
% (Aguilar-Rodríguez & Wagner, GBE, 2018)
blackList = {
    'cpd00067',... % H+
    'cpd00001',... % H2O
    'cpd00002',... % ATP
    'cpd00009',... % Phosphate
    'cpd00008',... % ADP
    'cpd00006',... % NADP
    'cpd00005',... % NADPH
    'cpd00003',... % NAD
    'cpd00004',... % NADH
    'cpd00012',... % PPi
    'cpd00011',... % CO2
    'cpd00010',... % CoA
    'cpd11493',... % ACP --> not in Aguilar-Rodríguez & Wagner, GBE, 2018
    'cpd00018',... % AMP
    'cpd00007',... % O2 --> not in Aguilar-Rodríguez & Wagner, GBE, 2018
    'cpd00013',... % NH3 --> not in Aguilar-Rodríguez & Wagner, GBE, 2018
    };

blackList = translateIDs(blackList, 'met',[], 'ModelSEED', 'MNXref');

% Cofactors -> Coenzymes
% (KEGG br08001: https://www.genome.jp/kegg-bin/get_htext#C42)
coFactorsKEGG = {
    'C00068';... %  TPP
    'C00016';... %  FAD
    'C00061';... %  FMN
    'C00003';... %  NAD
    'C00006';... %  NADP
    'C00010';... %  CoA
    'C00018';... %  PLP (Pyridoxal phosphate)
    'C00101';... %  THF (Tetrahydrofolate)
    'C00194';... %  Cobamamide
    'C00002';... %  ATP
    'C00063';... %  CTP
    'C00019';... %  SAM (S-Adenosyl-L-methionine)
    'C00029';... %  UDP-glucose
    'C04628';... %  Coenzyme B
    'C05777';... %  Coenzyme F430
    'C03576';... %  Coenzyme M
    'C11378';... %  Coenzyme Q10
    'C00051';... %  Glutathione
    'C00862';... %  Methanofuran
    'C00053';... %  PAPS (3-phosphoadenylylsulfate)
    'C00272';... %  Tetrahydrobiopterin
    'C01217';... %  5,6,7,8-Tetrahydromethanopterin
    'C15670';... %  Heme A
    'C00032';... %  Heme B
    'C15672';... %  Heme O
    'C00120';... %  Biotin
    'C00725';... %  Lipoate
    'C18237';... %  Molybdenum cofactor
    'C00113';... %  PQQ (Pyrroloquinoline-quinone)
    };

coFactors = translateIDs(coFactorsKEGG, 'met', [], 'KEGG', 'MNXref');

% Evaluate the draft models
analyzeRAVENModels
analyzeKBaseModels
analyzeCarveMeModel
analyzeAuReMeModels


%% translate gene IDs for KBase and AuReMe models
translation_of_gene_IDs

%% remove all exchange and biomass reactions from the models
remove_biomass_and_exchange_rxns

%% add metabolite formulae
add_formulae_to_models

%% adapt RAVEN GPR to COBRA standard
RAVEN_GPR_to_COBRA_format

%% Merge models from different approaches
merge_metabolic_models

%% add a universal prokaryotic biomass to all models
add_universal_biomass

%% Iterative gap filling
d = datestr(floor(now));
diary(['data/gap-filling/output-' d '.txt']);

% individual gap filling
gap_fill_individual_models

% conditional gap filling (COMMIT)
run_iterative_gap_filling
