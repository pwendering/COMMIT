% Impact of exchanged metabolites on biomass flux
options

% set COBRA solver to CPLEX
changeCobraSolver('ibm_cplex','all',0);

% medium that has been used for gap filling
load(mediumFile)

% taxonomic classification of At-SPHERE OTUs
% taxonomyFile = 'data/genomes/At-SPHERE-phyla.txt';
% taxonomyFile = 'data/genomes/At-SPHERE-classes.txt';
taxonomyFile = 'data/genomes/At-SPHERE-families.txt';

% load model workspace
subFolders = {'KBase', 'no_CarveMe', 'all'};

for E = experiments
    E = char(E);
    for F = subFolders
        
        F = char(F);
        
        fileBaseName = ['figures/exchanged_metabolites/reduction_biomass/',...
                habitat, '_', E, '_', F, '_biomass_'];
            
        fprintf('\n%s %s...\n', E, F)
        
        if ~exist([fileBaseName, 'export.txt'], 'file')
            
            modelWorkspace = ...
                fullfile('data/gap-filling/iterative',...
                    habitat, F, E);
            
            load(modelWorkspace, 'GF')
            n = numel(GF);
            
            %% extract exported metabolites from sink and uptake reactions that have been added
            % during iterative gap filling
            exported_per_model = cell(n, 1);
            imported_per_model = cell(n, 1);
            
            for i=1:n
                
                % find exported metabolites sink reactions
                exported = regexp(GF{i}.rxns(contains(GF{i}.rxns, 'sink_')),...
                    'MNXM\d+\[e]', 'match');
                exported_per_model{i} = [exported{:}]'; clear exported
                
                % find imported metabolites from exchange reactions
                imported = regexp(GF{i}.rxns(contains(GF{i}.rxns, 'EX_')),...
                    'MNXM\d+\[e]', 'match');
                imported_per_model{i} = [imported{:}]';
                
            end
            
            exported = unique(vertcat(exported_per_model{:}));
            imported = unique(vertcat(imported_per_model{:}));
            
            %% find optimal biomass flux values
            opt = zeros(n,1);
            for i=1:n
                v = optimizeCbModel(GF{i});
                opt(i) = v.f;
            end
            
            %% block uptake reactions for every metabolite separately
            bio_uptake = zeros(n, numel(imported));
            for i=1:numel(imported)
                rxnID = ['EX_', imported{i}];
                for j=1:n
                    % find the reaction index
                    idx = strcmp(GF{j}.rxns, rxnID);
                    % if the reaction is contained in the model, block it by setting
                    % upper and lower bounds to zero
                    if sum(idx) ~= 0
                        lb = GF{j}.lb(idx);
                        ub = GF{j}.ub(idx);
                        GF{j}.lb(idx) = 0;
                        GF{j}.ub(idx) = 0;
                        v = optimizeCbModel(GF{j});
                        bio_uptake(j,i) = v.f;
                        GF{j}.lb(idx) = lb;
                        GF{j}.ub(idx) = ub;
                    else
                        bio_uptake(j,i) = opt(j);
                    end
                    
                end
            end
            
            %% block sink reactions for every metabolite separately
            bio_sink = zeros(n, numel(exported));
            for i=1:numel(exported)
                rxnID = ['sink_', exported{i}];
                for j=1:n
                    % find the reaction index
                    idx = strcmp(GF{j}.rxns, rxnID);
                    % if the reaction is contained in the model, block it by setting
                    % upper and lower bounds to zero
                    if sum(idx) ~= 0
                        lb = GF{j}.lb(idx);
                        ub = GF{j}.ub(idx);
                        GF{j}.lb(idx) = 0;
                        GF{j}.ub(idx) = 0;
                        v = optimizeCbModel(GF{j});
                        bio_sink(j,i) = v.f;
                        GF{j}.lb(idx) = lb;
                        GF{j}.ub(idx) = ub;
                    else
                        bio_sink(j,i) = opt(j);
                    end
                    
                end
            end
            
            %% calculate the ratio to the optimal biomass
            ratio_sink = round(bsxfun(@rdivide, bio_sink, opt), 4);
            ratio_uptake = round(bsxfun(@rdivide, bio_uptake, opt), 4);
            
            model_ids = cellfun(@(x)strtok(x.id, '_'), GF, 'un', 0);
            
            %% write table to file
            writetable(array2table(ratio_sink, 'VariableNames', strtok(exported, '['),...
                'RowNames', model_ids),...
                [fileBaseName, 'export.txt'], 'Delimiter', '\t',...
                'WriteVariableNames', true, 'WriteRowNames', true)
            
            writetable(array2table(ratio_uptake, 'VariableNames', strtok(imported, '['),...
                'RowNames', model_ids),...
                [fileBaseName, 'import.txt'], 'Delimiter', '\t',...
                'WriteVariableNames', true, 'WriteRowNames', true)
        end
        
        %% classify metabolites using the KEGG BRITE hierarchy
        exported = readtable([fileBaseName, 'export.txt'],...
            'Delimiter', '\t', 'ReadVariableNames', true, 'ReadRowNames', true);
        exported = exported.Properties.VariableNames';
        
        imported = readtable([fileBaseName, 'import.txt'],...
            'Delimiter', '\t', 'ReadVariableNames', true, 'ReadRowNames', true);
        imported = imported.Properties.VariableNames';       
        
        brite_uptake = map2KEGGBrite(imported, briteFile);
        for i=1:numel(brite_uptake)
            if ~isempty([brite_uptake{i}{:}]) 
                brite_uptake(i) = brite_uptake{i}(1);
            else
                brite_uptake(i) = {'Other'};
            end
        end
        
        brite_export = map2KEGGBrite(exported, briteFile);
        for i=1:numel(brite_export)
            if ~isempty([brite_export{i}{:}])
                brite_export(i) = brite_export{i}(1);
            else
                brite_export(i) = {'Other'};
            end
        end
        brite_export(cellfun(@isempty,brite_export)) = {'Other'};
        
        writetable(cell2table([strtok(exported, '['), brite_export],...
            'VariableNames', {'ID', 'brite'}),...
            [fileBaseName, 'export_brite.txt'],...
            'WriteVariableNames', true, 'Delimiter', '\t')
        
        writetable(cell2table([strtok(imported, '['), brite_uptake],...
            'VariableNames', {'ID', 'brite'}),...
            [fileBaseName, 'import_brite.txt'],...
            'WriteVariableNames', true, 'Delimiter', '\t')
        
        
    end
end
