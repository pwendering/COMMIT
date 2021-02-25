% Run all evaluations

fprintf('##################################################################################\n')
fprintf('############### Starting Evaluation of genome-scale reconstructions ##############\n')
fprintf('##################################################################################\n\n\n')

tic
analyzeKBaseModels
analyzeRAVENModels
analyzeCarveMeModels
toc