% parse the MNXref universal reaction database and add metabolite
% permeability information to construct a gap-filling database

%% Preparation
% load options
options

compartments = {'c', 'e'};

blackList = {'BIOMASS',...
    'MNXM517',... % Light/hnu/photon
    'MNXM639',... % peptide
    };

% MNXref universal database filtered for balanced reactions
databaseFile = fullfile(topDir, 'data/gap-filling/database/reaction_MNXref_balanced.lst');

% parse database file
dbModel_MNXref_balanced = prepareFastGapFilling(databaseFile, compartments, blackList);

% save the output as a matlab workspace
save(topDir, 'data/gap-filling/database/Universal-model-MNXref-balanced.mat', 'dbModel_MNXref_balanced')
add_names_and_permeability_to_db




