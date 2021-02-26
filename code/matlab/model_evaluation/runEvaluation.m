% Run all evaluations
options; clearvars -except topDir

fprintf('##################################################################################\n')
fprintf('############### Starting Evaluation of genome-scale reconstructions ##############\n')
fprintf('##################################################################################\n\n\n')

tic
analyzeKBaseModels
analyzeRAVENModels
analyzeCarveMeModels
analyzeAuReMeModels
toc