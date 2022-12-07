% save the COBRA models as smbl
options

% set COBRA solver to CPLEX
changeCobraSolver('ibm_cplex','all',0);

experiment = 'Schlaeppi';

modelFile = fullfile('data/models/consensus/Soil_consensus_models_biomass.mat');
tmp = load(modelFile);
merged_models = tmp.merged_models;
clear tmp

% output directory
outDir = fullfile('data/models/consensus', habitat, '_sbml');

% create directory if it does not exist
if ~exist(outDir, 'dir')
    mkdir(outDir);
end

% processs model files
m_files = dir(outDir);
m_files = regexp({m_files.name}, '.*\.xml\.gz$', 'match');
m_files = [m_files{:}]';
m_files = strtok(m_files, '.');
n_blocked = zeros(numel(merged_models),1);

% fix the model fields for writing
for i=1:numel(merged_models)
    if ~ismember(merged_models{i}.id, m_files)
        tic
        disp(merged_models{i}.id)
        model = merged_models{i};
        % fix different types in fields 'rxnNotes', 'rxnNames' and 'EC'
        for j=1:numel(model.rxnNotes)
            if ~iscellstr(model.rxnNotes(j))
                model.rxnNotes(j) = cellstr(strjoin(model.rxnNotes{j}, ';'));
            end
            
            if ~iscellstr(model.rxnNames(j))
                try
                    model.rxnNames(j) = cellstr(strjoin(model.rxnNames{j}, ';'));
                catch
                    % biomass reaction name is vertical char array instead of
                    % horizontal
                    model.rxnNames(j) = cellstr(model.rxnNames{j}{1}');
                end
            end
            model.rxns = strrep(model.rxns,' ','_');
            
            if ~iscellstr(model.EC(j))
                model.EC(j) = cellstr(strjoin(model.EC{j}, ';'));
                model.EC(j) = regexprep(model.EC(j), '\|', ';');
            end
            
        end
        
        % make reaction IDs unique
        for j=1:numel(model.rxns)
            rxn = model.rxns(j);
            count = 1;
            while sum(ismember(model.rxns, rxn)) > 1
                count = count + 1;
                rxn = regexprep(rxn, '_\d$', '');
                rxn = cellstr(strcat(rxn{:}, '_',  num2str(count)));
                model.rxns(j) = rxn;
            end
        end
        
        
        % fix different types in fields 'metNames'
        for j=1:numel(model.metNames)
            if ~iscellstr(model.metNames(j))
                model.metNames(j) = model.metNames{j};
            end
        end
        model.mets = strrep(model.mets,' ','_');
        
        % format metabolite formulae
        for j=1:numel(model.metFormulas)
            % remove charge symbols from chemical formulae
            model.metFormulas(j) = regexprep(model.metFormulas(j), '[\-\+]+', '');
            % remove (.*)n from generic metabolites
            model.metFormulas(j) = regexprep(model.metFormulas(j), '\(.+\)[nx]+', '');
        end
                
        % format gene IDs
        model.genes = regexprep(model.genes, '^.*\|', '');
        
        % change EC field name to COBRA standard
        model.rxnECNumbers = model.EC;
        model = rmfield(model, 'EC');
        
        % adapt model descriptor field to COBRA standard
        model.modelID = model.id;
        model.modelName = model.description;
        model = rmfield(model, {'id', 'description'});
        
        % convert the model back to reversible
        model = convertModelToReversible(model);
        
        % NGAM reaction
        model = addReaction(model, 'MNXR151913_c',...
            'reactionName', 'non-growth associated maintenance reaction',...
            'reactionFormula', 'MNXM2[c] + 1 MNXM3[c] -> MNXM1[c] + MNXM7[c] + MNXM9[c]',...
            'printLevel', 0);
        
        model.metNames(findMetIDs(model, 'BIOMASS[c]')) = {'biomass'};
        model.rxnNames(findRxnIDs(model, 'BIOMASS_Reaction')) = {'universal_biomass_reaction'};
        model.rxns(findRxnIDs(model, 'BIOMASS_Reaction')) = {'growth_c'};
        
        % format model ids, add annotation and osensestr
        model.rxnMetaNetXID = strtok(model.rxns, '_');
        model.rxnMetaNetXID(cellfun(@isempty,...
            regexp(model.rxnMetaNetXID, 'MNXR\d+'))) = {''};
        comps = regexp(model.rxns, '_[ec][_0-9r]*$', 'match');
        idx = ~cellfun(@isempty, comps);
        comps = [comps{:}]';
        comps = regexprep(comps, '[_0-9r]+$', '');
        model.rxns = cellstr(strcat('r', num2str((1:numel(model.rxns))', '%04.f')));
        model.rxns(idx) = strcat(model.rxns(idx), comps);
        model.rxns(~idx) = strcat(model.rxns(~idx), '_e');
        
        model.metMetaNetXID = strtok(model.mets, '[');
        model.metMetaNetXID(cellfun(@isempty,...
            regexp(model.metMetaNetXID, 'MNXM\d+'))) = {''};
        comps = regexp(model.mets, '\[[ec]\]$', 'match');
        model.mets = cellstr(strcat('m', num2str((1:numel(model.mets))', '%04.f')));
        model.mets = strtok(model.mets, comps);
        model.osenseStr = 'max';
        
        % remove field subSystems
        model = rmfield(model, 'subSystems');
        
        % remove grRules field
        model = rmfield(model, 'grRules');
        
        model = orderModelFields(model);
        
        % write and compress sbml file
        disp('writing sbml file')
        writeCbModel(model, 'format', 'sbml', 'filename',fullfile(outDir, model.modelID));
        % gzip(fullfile(outDir, [model.modelID, '.xml']));
        fprintf('Finished! Time: %.2f min\n\n', toc/60)
        
        merged_models{i} = model;
    end
end