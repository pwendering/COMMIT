%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare sizes of gap-filling solutions per experiment for
% - iterative run with CarveMe models
% - iterative run without CarveMe models
% - iterative run with only KBase models
% - individual gap fillings with CarvMe models
% - individual gap fillings without CarveMe models
% - individual gap fillings with only KBase models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

habitat = 'Soil';

experiments = {'Schlaeppi'};

iterativeGfDir = fullfile('/stud/wendering/Masterthesis/DATA/Gap-filling/iterative',...
    habitat);
consensusDir = '/stud/wendering/Masterthesis/DATA/Consensus_models';

% directory to save the table for figure
figrDir = '/stud/wendering/Masterthesis/FIGURES/added_reactions';

n_comp = 6;

for i=1:numel(experiments)
    
    load(fullfile(consensusDir,...
        strcat(habitat, '_consensus_models_biomass_gf')), 'GF')
    added_reactions = cell(n_comp, numel(GF));
    added_reactions_no_exc = cell(n_comp, numel(GF));
    row_names = {'KBase_iterative', 'no_CarveMe_iterative', 'all_iterative',...
        'KBase_individual', 'noCarveMe_individual', 'all_individual'};
    model_ids = cellfun(@(x)strtok(x.id, '_'), GF,...
        'UniformOutput', false);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    row = 1;
    % load results from iterative gap filling with only KBase models
    load(fullfile(iterativeGfDir, 'KBase', experiments{i}), 'GF')
    % find added reactions
    n = numel(GF);
    for j=1:n
        
        added = GF{j}.rxns(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
        idx = strcmp(strtok(GF{j}.id, '_'), model_ids);
        added_reactions{row, idx} = added;
        
        added = added(~contains(added, 'EX_'));
        added_reactions_no_exc{row, idx} = added;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    row = 2;
    % load results from iterative gap filling without CarveMe models
    load(fullfile(iterativeGfDir,  'no_CarveMe', [experiments{i}]), 'GF')
    % find added reactions
    n = numel(GF);
    for j=1:n
        added = GF{j}.rxns(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
        idx = strcmp(strtok(GF{j}.id, '_'), model_ids);
        added_reactions{row, idx} = added;
        
        added = added(~contains(added, 'EX_'));
        added_reactions_no_exc{row, idx} = added;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    row = 3;
    % load results from iterative gap filling with CarveMe models
    load(fullfile(iterativeGfDir, 'all', experiments{i}), 'GF')
    % find added reactions
    n = numel(GF);
    for j=1:n
        added = GF{j}.rxns(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
        idx = strcmp(strtok(GF{j}.id, '_'), model_ids);
        added_reactions{row, idx} = added;
        
        added = added(~contains(added, 'EX_'));
        added_reactions_no_exc{row, idx} = added;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % OTU IDs contained in current community
    otu_ids = cellfun(@(x)strtok(x.id, '_'), GF,...
        'UniformOutput', false);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    row = 4;
    % find added reactions in individually-gap-filled models with only
    % KBase draft models
    load(fullfile('/stud/wendering/Masterthesis/DATA/models_KBase',habitat,...
        strcat(habitat, '_models_biomass_gf')), 'GF')
    n = numel(GF);
    for j=1:n
        added = GF{j}.rxns(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
        idx = strcmp(strtok(GF{j}.id, '_'), model_ids);
        added_reactions{row, idx} = added;
        
        added = added(~contains(added, 'EX_'));
        added_reactions_no_exc{row, idx} = added;
    end; clear added idx n
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    row = 5;
    % find added reactions in individually-gap-filled models without
    % CarveMe models
    load(fullfile(consensusDir,...
        strcat(habitat, '_consensus_models_noCarveMe_biomass_gf')), 'GF')
    n = numel(GF);
    for j=1:n
        added = GF{j}.rxns(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
        idx = strcmp(strtok(GF{j}.id, '_'), model_ids);
        added_reactions{row, idx} = added;
        
        added = added(~contains(added, 'EX_'));
        added_reactions_no_exc{row, idx} = added;
    end; clear added idx n
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    row = 6;
    % find added reactions in individually-gap-filled models with CarveMe
    % models
    load(fullfile(consensusDir,...
        strcat(habitat, '_consensus_models_biomass_gf')), 'GF')
    n = numel(GF);
    for j=1:n
        added = GF{j}.rxns(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
        idx = strcmp(strtok(GF{j}.id, '_'), model_ids);
        added_reactions{row, idx} = added;
        
        added = added(~contains(added, 'EX_'));
        added_reactions_no_exc{row, idx} = added;
    end; clear added idx n
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % take the subset of OTUs found in the experiment
    idx = ismember(model_ids, otu_ids);
    added_reactions(:, ~idx) = [];
    added_reactions_no_exc(:, ~idx) = [];
    
    n_added = cellfun(@numel, added_reactions);
    n_added_no_exc = cellfun(@numel, added_reactions_no_exc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save results
    writetable(array2table(n_added,...
        'VariableNames', otu_ids, 'RowNames', row_names),...
        fullfile(figrDir, [habitat, '_', experiments{i}, '_gf_sol_size.txt']),...
        'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');
    
    writetable(array2table(n_added_no_exc,...
        'VariableNames', otu_ids, 'RowNames', row_names),...
        fullfile(figrDir, [habitat, '_', experiments{i}, '_gf_sol_size_no_exc.txt']),...
        'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');
    
    % calculate the Jaccard index between each pair of gap-filling
    % solutions for every model
    
    dist_mat_array = cell(size(added_reactions, 2), 1);
    dist_mat_array_no_exc = cell(size(added_reactions, 2), 1);
    % loop over models
    for j=1:size(added_reactions, 2)
        tmp = double(size(added_reactions, 1));
        tmp_no_exc = double(size(added_reactions, 1));
        % loop over methods
        for k=1:size(added_reactions, 1)
            % loop over methods
            for l=1:size(added_reactions, 1)
                tmp(k, l) = 1 - numel(intersect(vertcat(added_reactions{k, j}),...
                    vertcat(added_reactions{l, j}))) / ...
                    numel(union(vertcat(added_reactions{k, j}),...
                    vertcat(added_reactions{l, j})));
                tmp_no_exc(k, l) = 1 - numel(intersect(vertcat(added_reactions_no_exc{k, j}),...
                    vertcat(added_reactions_no_exc{l, j}))) / ...
                    numel(union(vertcat(added_reactions_no_exc{k, j}),...
                    vertcat(added_reactions_no_exc{l, j})));
            end
        end
        dist_mat_array{j} = tmp;
        dist_mat_array_no_exc{j} = tmp_no_exc;
    end
    for j=1:numel(dist_mat_array)
        dist_mat_array{j}(isnan(dist_mat_array{j})) = 1;
        dist_mat_array_no_exc{j}(isnan(dist_mat_array_no_exc{j})) = 1;
    end
    dist_added_reactions = STATIS(dist_mat_array);
    dist_added_reactions_no_exc = STATIS(dist_mat_array_no_exc);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % save results
    writetable(array2table(dist_added_reactions,...
        'VariableNames', row_names, 'RowNames', row_names),...
        fullfile(figrDir, [habitat, '_', experiments{i}, '_dist_methods.txt']),...
        'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');
    
    writetable(array2table(dist_added_reactions_no_exc,...
        'VariableNames', row_names, 'RowNames', row_names),...
        fullfile(figrDir, [habitat, '_', experiments{i}, '_dist_methods_no_exc_rxns.txt']),...
        'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');
end