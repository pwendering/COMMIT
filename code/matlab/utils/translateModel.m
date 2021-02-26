function translatedModel = translateModel(model, source, target, translationDB, complete, verbose)
%% translatedModel = translateModel(model, source, target, translationDB, complete, verbose)
% Translates a model into a different name space. The underlying table for
% translation is either found in the working directory or in a directory on
% the MATLAB path.
% Input:                    
%           struct model:               metabolic model containing the
%                                       fields 'rxns' and 'mets'
%           char source:                source namespace (see translateIDs
%                                       function for info)
%           char target:                target namespace
%           struct translationDB:       contains fields 'metTab', 'rxnTab'
%           logical complete:           adds all available identifiers to
%                                       the model if 1, default: 0
%           logical verbose (optional): if true, print warnings and
%                                       progress statements (default: true)
% Output:   struct translatedModel:     model with translated field and
%                                       additional fields if complete==1
%           cellstr unmappedRxns:       array that contains the ids of
%                                       reactions of the source namespace 
%                                       that could not be mapped
%           cellstr unmappedMets:       unmapped metabolite source ids

translatedModel = model;
if nargin < 4
    error('No translation table given')
elseif nargin < 5
    complete = 0;
    verbose = 1;
elseif nargin < 6 || ~islogical(verbose)
    verbose = 1;
elseif isempty(target) || ~ischar(source)
    warning('Parameter complete has been set to true since you have not given a valid source namespace.')
    complete = 1;    
elseif isempty(source)
    error('You have not given a valid input namespace.')
elseif ~ismember(complete, [0 1 true false])
     warning('Parameter complete has to be in [0,1,true,false]')
     complete = 0;
end

% Extract the reaction ids, suffixes and compartment ids
rxn_comps = regexp(model.rxns, '_', 'split');
rxn_comps = cellfun(@(x)x{end}, rxn_comps, 'UniformOutput', false);
rxns = regexprep(model.rxns, '_.$', '');
% get the ids of biomass reaction ans boundary reactions
ex_rxns = findExcRxns(model);
biomass = logical(model.c);

% Metabolite identifiers also contain compartment information '[x]'
mets = strtok(model.mets, '[');
met_comps = regexp(model.mets, '\[.\]', 'match');
met_comps = vertcat(met_comps{:});

if complete
    if isempty(target)
        return
        % dont translate, only add fields
    else
        return
    % translate to all options
    end
else
    if verbose
        fprintf('\nReactions: %s ==> %s\n', source, target)
    end
    
    % translate the reaction iss to the target namespace
    new_rxns = translateIDs(rxns, 'rxn', translationDB.rxnTab, source, target, verbose);
        
    % add prefixes compartment identifiers to the reaction ids (except for the
    % biomass reaction and boundary reactions)
    idx_nm = cellfun('isempty', new_rxns);
    new_rxns = strcat(new_rxns, '_', rxn_comps);
    
    % replace the non-translated ids with the original ids
    
    idx_nm = logical(ex_rxns + biomass + idx_nm);
    new_rxns(idx_nm) = model.rxns(idx_nm);
    
    if verbose
        fprintf('\nMetabolites: %s ==> %s\n', source, target)
    end
    
    % translate the metabolite ids to the target namespace
    new_mets = translateIDs(mets, 'met', translationDB.metTab,  source, target, verbose);
    % replace the non-translated ids with the original ids
    idx_nm = cellfun('isempty', new_mets);
    new_mets(idx_nm) = mets(idx_nm);
    % add compartment idetifiers to the metabolite ids
    new_mets = strcat(new_mets, met_comps);
    
    translatedModel.rxns = new_rxns;
    translatedModel.mets = new_mets;
    
end
    
end