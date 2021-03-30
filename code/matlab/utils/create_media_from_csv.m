options

%% Leaf gap-fill medium
% medium = readtable('data/media/Leaf/Leaf-media-combined-augmented.csv');
% medium = translateIDs(medium.compounds, 'met', [], 'ModelSEED', 'MNXref');
% medium = strcat(medium, '[e]');
% save('data/media/Leaf/leaf-medium.mat', 'medium')

%% Root gap-fill medium
% medium = readtable('data/media/Root/Root-media-combined-augmented.csv');
% medium = translateIDs(medium.compounds, 'met', [], 'ModelSEED', 'MNXref');
% medium = strcat(medium, '[e]');
% save('data/media/Root/root-medium.mat', 'medium')

%% Soil gap-fill medium
% medium = readtable('data/media/Soil/auxo-media/combined-media.csv');
% medium = translateIDs(medium.Var2, 'met', [], 'ModelSEED', 'MNXref');
% medium = strcat(medium, '[e]');
% save('data/media/Soil/soil-medium.mat', 'medium')
% 
% medium = readtable('data/media/Soil/soil-medium-lab.csv',...
%     'ReadVariableNames', false);
% medium = translateIDs(medium.Var1, 'met', [], 'ModelSEED', 'MNXref');
% medium = strcat(medium, '[e]');
% save('data/media/Soil/soil-medium-lab.mat', 'medium')

%% minimal medium
medium = readtable('data/media/minimal-medium.csv',...
    'ReadVariableNames', false);
medium = translateIDs(medium.Var1, 'met', [], 'ModelSEED', 'MNXref');
medium = strcat(medium, '[e]');
save('data/media/minimal-medium.mat', 'medium')