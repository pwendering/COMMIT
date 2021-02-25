% Generating input files for COMMGEN
load('/stud/wendering/Masterthesis/DATA/models_KBase/Soil/Soil_models_metFormulas.mat')
model = models{end}; clear models
model_KB = model;
model = removeRxns(model, model.rxns(contains(model.rxns, 'EX_')));
outDir = '/stud/wendering/Masterthesis/DATA/COMMGEN/COMMGEN_test/Soil811_KB';
generateCOMMGENInputFiles(model, outDir)


load('/stud/wendering/Masterthesis/DATA/models_CarveMe/Soil/Soil_models_metFormulas.mat')
model = models{end}; clear models
model_CM = model;
model = removeRxns(model, model.rxns(contains(model.rxns, 'EX_')));
outDir = '/stud/wendering/Masterthesis/DATA/COMMGEN/COMMGEN_test/Soil811_CM';
generateCOMMGENInputFiles(model, outDir)

clear model

COMMGEN({'/stud/wendering/Masterthesis/DATA/COMMGEN/COMMGEN_test/Soil811_KB',...
    '/stud/wendering/Masterthesis/DATA/COMMGEN/COMMGEN_test/Soil811_CM'})


load('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat')

merged = mergeModels({model_KB, model_CM}, dbModel_MNXref_balanced);