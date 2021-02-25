
%% Leaf gap-fill medium
medium = readtable('/stud/wendering/Masterthesis/DATA/media/Leaf/Leaf-media-combined-augmented.csv');
medium = translateIDs(medium.compounds, 'met', [], 'ModelSEED', 'MNXref');
medium = strcat(medium, '[e]');
save('/stud/wendering/Masterthesis/DATA/media/Leaf/leaf-medium.mat', 'medium')

%% Root gap-fill medium
medium = readtable('/stud/wendering/Masterthesis/DATA/media/Root/Root-media-combined-augmented.csv');
medium = translateIDs(medium.compounds, 'met', [], 'ModelSEED', 'MNXref');
medium = strcat(medium, '[e]');
save('/stud/wendering/Masterthesis/DATA/media/Root/root-medium.mat', 'medium')

%% Soil gap-fill medium
medium = readtable('/stud/wendering/Masterthesis/DATA/media/Soil/auxo-media/combined-media.csv');
medium = translateIDs(medium.Var2, 'met', [], 'ModelSEED', 'MNXref');
medium = strcat(medium, '[e]');
save('/stud/wendering/Masterthesis/DATA/media/Soil/soil-medium.mat', 'medium')

medium = readtable('/stud/wendering/Masterthesis/DATA/media/Soil/soil-medium-lab.csv',...
    'ReadVariableNames', false);
medium = translateIDs(medium.Var1, 'met', [], 'ModelSEED', 'MNXref');
medium = strcat(medium, '[e]');
save('/stud/wendering/Masterthesis/DATA/media/Soil/soil-medium-lab.mat', 'medium')

%% minimal medium
medium = readtable('/stud/wendering/Masterthesis/DATA/media/minimal-medium.csv',...
    'ReadVariableNames', false);
medium = translateIDs(medium.Var1, 'met', [], 'ModelSEED', 'MNXref');
medium = strcat(medium, '[e]');
save('/stud/wendering/Masterthesis/DATA/media/minimal-medium.mat', 'medium')