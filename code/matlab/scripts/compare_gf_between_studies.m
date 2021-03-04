% compare the added reactions between gap fillings
options

iterativeGfDir = fullfile('data/gap-filling/iterative', habitat);

singleGfFile = fullfile('data/models/kbase', habitat, [habitat, '_models_biomass_gf.mat']);
% singleGfFile = fullfile('data/models/consensus', [habitat, '_consensus_models_biomass_gf']);
% singleGfFile = fullfile('data/models/consensus', [habitat, '_consensus_models_noCarveMe_biomass_gf']);

% medium that has been used for gap filling
load(mediumFile)

% directory to save the table for figure
figOutDir = 'figures';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Compare the IDs of added reactions and exchanged metabolites between studies
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% load results of individual gap fillings
load(singleGfFile, 'GF')


model_ids = cellfun(@(x)strtok(x.id, '_'), GF,...
    'UniformOutput', false);
added_reactions = cell(numel(experiments), numel(GF));

imported = cell(numel(experiments), numel(GF));
exported = cell(numel(experiments), numel(GF));

for i=1:numel(experiments)
    
    % load results from iterative gap filling
    load(fullfile(iterativeGfDir, sub_dir, experiments{i}), 'GF')
    
    n = numel(GF);
    
    for j=1:n
        
        idx = strcmp(strtok(GF{j}.id, '_'), model_ids);
        
        % find added reactions
        added = GF{j}.rxns(cellfun(@(x)isequal(x,'gf'), GF{j}.rxnNotes));
        added_reactions{i, idx} = added;
                
        % find exported metabolites sink reactions
        exp = regexp(GF{j}.rxns(contains(GF{j}.rxns, 'sink_')),...
            'MNXM\d+\[e]', 'match');
        exp = [exp{:}]';
        exported{i, idx} = exp;
        
        % find imported metabolites from exchange reactions
        imp = regexp(GF{j}.rxns(contains(GF{j}.rxns, 'EX_')),...
            'MNXM\d+\[e]', 'match');
        imp = [imp{:}]';
        imp = setdiff(imp, medium);
        imported{i, idx} = imp;
        
    end; clear added idx n
    
    
    
end

% take the intersect of OTUs found in both the experiments
idx = all(~cellfun('isempty', added_reactions), 1);
added_reactions(:, ~idx) = [];
exported(:, ~idx) = [];
imported(:, ~idx) = [];
model_ids(~idx) = [];

% find the intersection and set differences between the studies for each
% model (for reactions, exported and imported metabolties)
shared_rxn = zeros(numel(experiments)+1, numel(model_ids));
shared_exp = zeros(numel(experiments)+1, numel(model_ids));
shared_imp = zeros(numel(experiments)+1, numel(model_ids));

for i=1:numel(model_ids)
    for j=1:numel(experiments)
        
        % reactions/metabolites of the current model and experiment
        idx_rest = true(numel(experiments), 1);
        idx_rest(j) = 0;
        
        R_current = added_reactions{j,i};
        R_rest = vertcat(added_reactions{idx_rest,i});
        
        E_current = exported{j,i};
        E_rest = vertcat(exported{idx_rest,i});
        
        I_current = imported{j,i};
        I_rest = vertcat(imported{idx_rest,i});
        
        % find the number of reactions/metabolites that is unique to the
        % current model+experiment
        shared_rxn(j, i) = numel(R_current) - numel(intersect(R_current, R_rest));
        shared_exp(j, i) = numel(E_current) - numel(intersect(E_current, E_rest));
        shared_imp(j, i) = numel(I_current) - numel(intersect(I_current, I_rest));
        
        % intersection of all experiments per model (last row)
        if j==1
            R_i = intersect(R_current, R_rest);
            E_i = intersect(E_current, E_rest);
            I_i = intersect(I_current, I_rest);
        else
            R_i = intersect(R_i, intersect(R_current, R_rest));
            E_i = intersect(E_i, intersect(E_current, E_rest));
            I_i = intersect(I_i, intersect(I_current, I_rest));
        end
    end
    
    shared_rxn(end, i) = numel(R_i);
    shared_exp(end, i) = numel(E_i);
    shared_imp(end, i) = numel(I_i);
end; clear R_current R_rest E_current E_rest I_current I_rest

% save results

tmpOutDir = fullfile(figOutDir, 'added_reactions');

% Reactions
writetable(array2table(shared_rxn, 'VariableNames', model_ids,...
    'RowNames', [experiments, {'Intersection'}]),...
    fullfile(tmpOutDir, [sub_dir,'_', habitat, '_comp_gf_exp.txt']),...
    'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');

tmpOutDir = fullfile(figOutDir, 'exchanged_metabolites');

% Exported metabolites
writetable(array2table(shared_exp, 'VariableNames', model_ids,...
    'RowNames', [experiments, {'Intersection'}]),...
    fullfile(tmpOutDir, [sub_dir,'_',habitat, '_exported.txt']),...
    'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');

% Imported metabolites
writetable(array2table(shared_imp, 'VariableNames', model_ids,...
    'RowNames', [experiments, {'Intersection'}]),...
    fullfile(tmpOutDir, [sub_dir,'_',habitat, '_imported.txt']),...
    'WriteVariableNames', true, 'WriteRowNames', true, 'Delimiter', '\t');