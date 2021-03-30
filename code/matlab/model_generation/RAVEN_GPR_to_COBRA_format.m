% convert RAVEN GPR rules to COBRA standard
options; clear

habitats = {'Soil', 'Leaf', 'Root'};

disp('-------------------------------------------------------------------')
disp('START')
disp('-------------------------------------------------------------------')

modelDir = 'data/models/raven/HMMer10E-50';

for i=1:numel(habitats)
    disp(habitats{i})
    disp('------')
    workspace = fullfile(modelDir, strcat(habitats{i}, '_models_metFormulas.mat'));
    load(workspace)
    for j=1:numel(models)
        model = models{j};
        fprintf('converting %s ...\n', models{j}.id)
        % replace every gene name by its index
        for k=1:numel(model.genes)
            id = model.genes{k};
            id = regexptranslate('escape', id);
            id = strcat(id, '\>');
            repl = strcat('x(', num2str(k), ')');
            model.rules = cellfun(@(x)regexprep(x, id, repl),...
                model.rules, 'UniformOutput', false);
        end
        model.rules =  cellfun(@(x)regexprep(x, 'or', '\|'),...
            model.rules, 'UniformOutput', false);
        models{j} = model;
        disp('done')
        disp('------')
    end
    workspace = fullfile(modelDir, strcat(habitats{i}, '_models_COBRA_GPR.mat'));
    save(workspace, 'models')
    clear models
    disp('-------------------------------------------------------------------')
end
