function consensusModel = mergeModels(models, database)
% Merge a list of metabolic models based on their reaction and metabolite
% identifiers and stoichiometric coefficients. All models have to be in the
% MetaNetX namespace.
%
% Input:
%           cell models:            array of structs (metabolic models)
%           logical dbCheck:        whether each reaction should be checked
%                                   against a reaction with the identical ID
%                                   in the provided database model
%           struct database:        database model to be checked against
%                                   (ensure mass balanced reactions!)
% Output:
%           struct consensusModel:  merged metabolic model
%
%

if nargin < 2
    error('You did not give enough arguments.')
end

models = reshape(models, numel(models), 1);
% Create an empty struct
consensusModel = models{1};
comps = models{1}.comps;
compNames = models{1}.compNames;

for i=2:numel(models)
    tic
    fprintf('Integrating model %d...', i)
    model = models{i};
    %% Add the model to the consensus model by extending every contained field
    fields_current = fieldnames(model);
    fields_consensus = fieldnames(consensusModel);
    
    for j=1:numel(fields_current)
        field = fields_current{j};
        % if field is present in the merged model, extend the field
        if field == 'S'
            consensusModel.S(end+1:end+size(model.S,1), end+1:end+size(model.S, 2)) = model.S;
        elseif isfield(consensusModel, field) && ~ismember(field, {'rules', 'genes'})
            try
                consensusModel.(field) = vertcat(consensusModel.(field), model.(field));
            catch
                consensusModel.(field) = '';
            end
        end
    end
    
    % remove all fields from the consensus model that are not contained in
    % the current model
    consensusModel = rmfield(consensusModel, setdiff(fields_consensus, fields_current));
    
    %% update gene rules
    for j=1:numel(models{i}.rxns)
        % indices of the associated genes
        gene_idx = cellfun(@str2double,regexp(models{i}.rules{j}, '[0-9]*', 'match'));
        % identifiers of these genes
        genes_tmp = models{i}.genes(gene_idx);
        if ~isempty(genes_tmp)
            % initialize the new GPR rule
            new_rule = {''};
            % extract all operators to be able to reconstruct the rule
            operators = regexpi(models{i}.rules{j}, '\||\&|and|or', 'match');
            for k=1:numel(genes_tmp)
                % for each gene, check whether it is already contained in
                % the merged model
                if ~ismember(genes_tmp(k), consensusModel.genes) && ~isempty(genes_tmp{k})
                    % if the gene is not contained, add the gene and the
                    % new index for the GPR rule is the number of genes
                    consensusModel.genes = vertcat(consensusModel.genes, genes_tmp(k));
                    new_idx = numel(consensusModel.genes);
                    % if the gene is contained, find the index
                elseif isempty(genes_tmp{k})
                    new_idx = 0;
                else
                    new_idx = find(strcmp(consensusModel.genes, genes_tmp{k}));
                end
                
                if new_idx~=0
                    
                    % check if the gene occurs multiple times in the
                    % consensus model
                    if numel(new_idx)>1
                        new_idx = new_idx(1);
                    end
                    new_idx = num2str(new_idx);
                    % build up the updated GPR rule
                    if isempty(new_rule{:})
                        % if first element, add it
                        new_rule = {strcat('x(', new_idx, ')')};
                    else
                        % if not first element, join the rules by the current
                        % operator between them
                        new_rule = {strjoin(vertcat(new_rule, strcat('x(', new_idx, ')')),...
                            char(strcat({' '}, operators{k-1}, {' '})))};
                    end
                end
            end
            % add the new rule to the merged model
            consensusModel.rules = vertcat(consensusModel.rules, new_rule);
        else
            % if no gene associated to the reaction, add an empty entry
            consensusModel.rules = vertcat(consensusModel.rules, {''});
        end
        clear gene_idx genes_tmp new_idx new_rule operators
    end
    
    %% remove duplicate metabolites
    consensusModel = removeDuplicateMets(consensusModel);
    
    %% remove duplicate reactions
    
    consensusModel = removeDuplicateRxns(consensusModel, database, 0);
    fprintf('(%.2fs)\n', toc)
end

%%  id and description

consensusModel.id = models{1}.id;
consensusModel.description = strcat(strtok(models{1}.description, '_'),...
    '_merged');
consensusModel.comps = comps;
consensusModel.compNames = compNames;

if isfield(consensusModel, 'grRules')
    consensusModel = rmfield(consensusModel, 'grRules');
end

end