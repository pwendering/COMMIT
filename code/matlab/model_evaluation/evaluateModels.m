function evaluateModels(models, recMethod, plotDir, resWorkspace)
%% evaluateModels(models, recMethod, plotDir, resWorkspace)
% Evaluation of models by calling the distance method functions and save
% the resulting matrices and figures
% Input:
%       cell models:            array of struct objects (metabolic models)
%       char recMethod:         method that was used to create the models 
%                               (will be used to name the figure files)
%       char plotDir:           path to the directory in which the figures
%                               should be saved
%       char res_Workspace:     name of the .mat workspace where the
%                               resulting matrices should be stored
% Output:
%       saves matrices in workspace and saves figures for each distance
%       measure

options; clear
% load tables required for some functions
tablesDir = 'data/tables';

ecTranslationTable = readtable(fullfile(tablesDir, 'corrected-EC-numbers.csv'),...
    'ReadVariableNames', false);

coFactorsTab = readtable(fullfile(tablesDir, 'cofactors_from_KEGG.csv'),...
    'Delimiter', '\t', 'ReadVariableNames', false);
coFactors = coFactorsTab.Var5; % MNXref
coFactorsNames = coFactorsTab.Var4;

if plotDir
    saveFigures = true;
else 
    saveFigures = false;
end

n = numel(models);

% convert the stoichiometric matrices to full
for j=1:n
    models{j}.S = full(models{j}.S);
end

%% Obtain the taxonomy names from the field 'description'
fprintf('\nProcessing taxonomy names...')
tax_names = cell(n,1);
tax_names(:) = {''};
for j=1:n
    name = models{j}.description;
    split = strsplit(name ,'_');
    count = 1;
    name = strcat(split{1},'_', num2str(count));
    while ismember(name, tax_names)
        count = count + 1;
        name = regexpi(name, '[a-z]*', 'match');
        name = strcat(name{:}, '_',  num2str(count));
    end
    tax_names{j} = name;
end
clear split name count
fprintf(' done.\n')

%% Re-arrange the model array according to species names
[tax_names, idx] = sort(tax_names);
models = models(idx);

%% Evaluation
fprintf('\nStarting the evaluation.\n\n')

tic

% SVD-based distance
fprintf('\t> SVD-based distance...')
SVD_matrix = svdDistance(models);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-SVD.png'));
plotDistance(SVD_matrix, tax_names, saveFigures, filename, 'SVD')
fprintf(' done. (%.0fs)\n', toc); tmp_toc = toc;

% Jaccard distance between the reaction sets
fprintf('\t> reaction Jaccard distance...')
rxnJaccard_matrix = rxnJaccardDist(models);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-rxnJaccard.png'));
plotDistance(rxnJaccard_matrix, tax_names, saveFigures, filename, 'reaction-Jaccard')
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

% Jaccard distance between the metabolite sets
fprintf('\t> metabolite Jaccard distance...')
metJaccard_matrix = metJaccardDist(models);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-metJaccard.png'));
plotDistance(metJaccard_matrix, tax_names, saveFigures, filename, 'metabolite-Jaccard')
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

% Difference in numbers of metabolic dead-ends
fprintf('\t> numbers of dead end metabolites...')
deadEndsNumber_matrix = deadEndsNumberDist(models);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-numberDeadEnds.png'));
plotDistance(deadEndsNumber_matrix, tax_names, saveFigures, filename, '#dead-ends')
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

% Jaccard distance between the set of dead-end metabolites
fprintf('\t> dead-end Jaccard distance...')
deadEndJaccard_matrix = deadEndsJaccardDist(models);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-deadEndsJaccard.png'));
plotDistance(deadEndJaccard_matrix, tax_names, saveFigures, filename, 'dead-ends-Jaccard')
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

% Jaccard distance between the EC numbers
fprintf('\t> E.C. number Jaccard distance...')
ecLevel = 4;
ecJaccard_matrix = ecJaccardDist(models, ecLevel);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-ecJaccard',num2str(ecLevel),'.png'));
plotDistance(ecJaccard_matrix, tax_names, saveFigures, filename, strcat('ec-Jaccard-lv',num2str(ecLevel)))
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

% Difference in abundance of EC numbers
fprintf('\t> abundances of E.C. numbers...')
ecAbundanceDist_matrix = ecAbundanceDist(models, ecTranslationTable.Var2);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-ecAbundance.png'));
plotDistance(ecAbundanceDist_matrix, tax_names, saveFigures, filename, 'abundance-of-EC-numbers')
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

% Distance based on usage of cofactors
fprintf('\t> co-factor usage...')
cfmat = cofactorUsage(models, coFactors);
fileName = fullfile(plotDir, strcat('heatmap-models-', recMethod, '-cofactor-abundance.png'));
f = figure('visible','off');
heatmap(coFactorsNames,tax_names,cfmat);
saveas(f, fileName)
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

fprintf('\t> co-factor abundance...')
coFactorDistance_matrix = coFactorDist(models, coFactors);
filename = fullfile(plotDir,strcat('heatmap-models-', recMethod, '-coFactors.png'));
plotDistance(coFactorDistance_matrix, tax_names, saveFigures, filename, 'Usage of cofactors')
fprintf(' done. (%.0fs)\n', toc-tmp_toc); tmp_toc = toc;

tev = toc;

% Save all matrices
matrix_objects = who('*matrix*');
save(resWorkspace, matrix_objects{:})
save(resWorkspace, 'tev', 'cfmat', '-append')
        
end