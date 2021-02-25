% Impact of exchanged metabolites on biomass flux

% initCobraToolbox(false)

habitat = 'Root';

experiments = {'Bulgarelli', 'Schlaeppi'};

% medium that has been used for gap filling
mediumFile = '/stud/wendering/Masterthesis/DATA/media/minimal-medium.mat';
load(mediumFile)

% taxonomic classification of At-SPHERE OTUs
% taxonomyFile = '/stud/wendering/Masterthesis/DATA/genomes/At-SPHERE-phyla.txt';
% taxonomyFile = '/stud/wendering/Masterthesis/DATA/genomes/At-SPHERE-classes.txt';
taxonomyFile = '/stud/wendering/Masterthesis/DATA/genomes/At-SPHERE-families.txt';

% load model workspace
subFolders = {'all'};%{'KBase', 'no_CarveMe', 'all'};

briteFile = '/net/calc1/srv/wendering/kegg/brite/briteHierarchy_ext.csv';
pathwayFile = '/net/calc1/srv/wendering/kegg/pathway/pathway_compounds_red.csv';

for E = experiments
    E = char(E);
    for F = subFolders
        
        F = char(F);
        
        fileBaseName = ['/stud/wendering/Masterthesis/FIGURES/exchanged_metabolites/reduction_biomass/',...
                habitat, '_', E, '_', F, '_biomass_'];
            
        fprintf('\n%s %s...\n', E, F)
        
        if ~exist([fileBaseName, 'export.txt'], 'file')
            
            modelWorkspace = ...
                fullfile('/stud/wendering/Masterthesis/DATA/Gap-filling/iterative',...
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
        
        %% classify metabolites
        exported = readtable([fileBaseName, 'export.txt'],...
            'Delimiter', '\t', 'ReadVariableNames', true, 'ReadRowNames', true);
        exported = exported.Properties.VariableNames';
        
        imported = readtable([fileBaseName, 'import.txt'],...
            'Delimiter', '\t', 'ReadVariableNames', true, 'ReadRowNames', true);
        imported = imported.Properties.VariableNames';
        
        % pathways
%         pathways_uptake = map2KEGGPathway(imported, pathwayFile);
%         pathways_uptake_all = vertcat(pathways_uptake{:});
%         pathways_uptake_unique = unique(pathways_uptake_all);
%         n_pathways_uptake = cellfun(@(x)sum(strcmp(pathways_uptake_all, x)),...
%             pathways_uptake_unique); clear pathways_uptake_all
%         %         [~,ia] = sort(n_pathways_uptake, 'descend');
%         for i=1:numel(pathways_uptake)
%             if numel(pathways_uptake{i}) > 1
%                 %                 pathways_uptake(i) = pathways_uptake{i}(1);
%                 pathways_uptake{i} = strjoin(pathways_uptake{i}, ';');
%             end
%         end
%         
%         pathways_export = map2KEGGPathway(exported, pathwayFile);
%         pathways_export_all = vertcat(pathways_export{:});
%         pathways_export_unique = unique(pathways_export_all);
%         n_pathways_export = cellfun(@(x)sum(strcmp(pathways_export_all, x)),...
%             pathways_export_unique); clear pathways_export_all
%         %         [~,ia] = sort(n_pathways_export, 'descend');
%         
%         for i=1:numel(pathways_export)
%             if numel(pathways_export{i}) > 1
%                 %                 pathways_export(i) = pathways_export{i}(1);
%                 pathways_export{i} = strjoin(pathways_export{i}, ';');
%             end
%         end
%         
%         writetable(cell2table([imported, pathways_uptake],...
%             'VariableNames', {'ID', 'pathway'}),...
%             [fileBaseName, 'import_pathway.txt'],...
%             'WriteVariableNames', true, 'Delimiter', '\t')
%         
%         writetable(cell2table([exported, pathways_export],...
%             'VariableNames', {'ID', 'pathway'}),...
%             [fileBaseName, 'export_pathway.txt'],...
%             'WriteVariableNames', true, 'Delimiter', '\t')
        
        
        % brite
        brite_uptake = map2KEGGBrite(imported, briteFile);
        for i=1:numel(brite_uptake)
            if numel(brite_uptake{i}) > 1
                brite_uptake(i) = brite_uptake{i}(1);
            end
        end
        
        brite_export = map2KEGGBrite(exported, briteFile);
        for i=1:numel(brite_export)
            if numel(brite_export{i}) > 1
                brite_export(i) = brite_export{i}(1);
            end
        end
        
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

%% Also: sum of all metabolites in the models together
