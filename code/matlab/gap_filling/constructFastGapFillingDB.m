%%%%%% runFastGapFilling

%% Preparation

compartments = {'c', 'e'};

blackList = {'BIOMASS',...
    'MNXM517',... % Light/hnu/photon
    'MNXM639',... % peptide
    };

% whole ModelSEED database
% databaseFile = '/net/calc1/srv/wendering/ModelSEED/reaction.lst';
% databaseFile = '/net/calc1/srv/wendering/mnxref/reaction_MNXref.lst';
databaseFile = '/net/calc1/srv/wendering/mnxref/reaction_MNXref_balanced.lst';
% databaseFile = '/net/calc1/srv/wendering/ModelSEED/reaction_KBase.lst';
% KBase gap filling database (ModelSEED extract)

% dbModel = prepareFastGapFilling(databaseFile, compartments, blackList);
% dbModel_KBase = prepareFastGapFilling(databaseFile, compartments, blackList);
% dbModel_MNXref = prepareFastGapFilling(databaseFile, compartments, blackList);


dbModel_MNXref_balanced = prepareFastGapFilling(databaseFile, compartments, blackList);
save('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat', 'dbModel_MNXref_balanced')
add_names_and_permeability_to_db;

% check 'MNXR124010'

% save('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref.mat', 'dbModel_MNXref')
% save('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-ModelSEED.mat', 'dbModel')
% save('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-ModelSEED.mat', 'dbModel_KBase', '-append')



