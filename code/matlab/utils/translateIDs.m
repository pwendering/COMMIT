function trList = translateIDs(idList, type, dbTable, source, target, verbose)
%% trList = translateIDs(idList, type, dbTable, source, target, verbose)
% Translate metabolite IDs from one namespace to another (ModelSEED, KEGG,
% MetaCyc, BiGG).
% Input:
%       cell idList:               array containing the ids to be translated
%       char type:                  either 'rxn' or 'met'
%       table dbTable:              table that contains the translation
%                                   of identifiers for the given type
%       char source:                source namespace ('ModelSEED',
%                                   'KEGG', 'MetaCyc', 'BiGG', 'MNXref',
%                                   (for rxns: 'Rhea', 'EC')
%                                   (for mets: 'ChEBI', 'NAMES')
%       char target:                target namespace
%       logical verbose (optional): if true, print warnings and
%                                   progress statements (default: true)
% Output:
%       cell trList:                array containing the translated IDs

if nargin < 6 || ~islogical(verbose)
    verbose = true;
end

options
clearvars -except topDir
dbDir = fullfile(topDir, 'data/tables/MNXref');

if ~iscellstr(idList)
    if ischar(idList)
        idList = {idList};
    elseif iscell(idList)
        idList = cellfun(@(x)char(x), idList, 'UniformOutput', false);
    else
        error('The input idList is either not given or not of type cellstr')
    end
end
    
if ~ischar(source) || ~ischar(target) || any(~ismember(type, ['rxn' 'met']))
    error('The source and/or target namespace or type definition is incorrect')
end

% Available namespaces
rxnSources = {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED', 'Rhea', 'EC'};
metSources = {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED', 'ChEBI', 'NAMES'};

%% Find matched using the MNXref database

% get the respective column that contains the desired names spaces
if type == 'met'
    sourceID = find(contains(metSources, source));
    targetID = find(contains(metSources, target));
else
    sourceID = find(contains(rxnSources, source));
    targetID = find(contains(rxnSources, target));
end

if isempty(sourceID) || isempty(targetID)
    error('The requested namespace is not available')
end

if ~isempty(dbTable) && numel(idList) <= 300
    trList = repmat({''}, numel(idList), 1);
    % Filter source  and target IDs by empty keys and keys that are definitely not
    % contained in the given list
    sourceIDs = dbTable.(source);
    targetIDs = dbTable.(target);
    sourceIDs = sourceIDs(contains(dbTable.(source), idList));
    targetIDs = targetIDs(contains(dbTable.(source), idList));
    clear dbTable;
    
    % ID at the beginning of the line followed by other IDs
    idList_suffix = strcat('^', idList, '\|');
    idx = cellfun(@(x)regexp(sourceIDs, x), idList_suffix, 'UniformOutput', false);
    for i=1:numel(idx)
        idx{i} = find(~cellfun('isempty', idx{i}));
    end
    clear idList_suffix
    
    % ID at the end of the line, preceded by other IDs
    empty = cellfun('isempty', idx);
    idList_prefix = strcat('\|', idList(empty), '$');
    idx_prefix = cellfun(@(x)regexp(sourceIDs, x), idList_prefix, 'UniformOutput', false);
    for i=1:numel(idx_prefix)
        idx_prefix{i} = find(~cellfun('isempty', idx_prefix{i}));
    end
    idx(empty) = idx_prefix; clear idx_prefix idList_prefix
    
    % ID in the middle of other IDs
    empty = cellfun('isempty', idx);
    idList_center = strcat('\|', idList(empty), '\|');
    idx_center = cellfun(@(x)regexp(sourceIDs, x), idList_center, 'UniformOutput', false);
    for i=1:numel(idx_center)
        idx_center{i} = find(~cellfun('isempty', idx_center{i}));
    end
    idx(empty) = idx_center; clear idx_center idList_center
    
    % field only contains the ID
    empty = cellfun('isempty', idx);
    idx(empty) = cellfun(@(x)find(strcmp(x, sourceIDs)), idList(empty), 'UniformOutput', false);
    
    for i=1:numel(idList)
        if ~isempty(idx{i}) && numel(idx{i}) == 1
            trList(i) = targetIDs(idx{i});
        else
            trList(i) = idList(i);
        end
    end
else
    dbFile = [dbDir, '/MNXref-',type,'-translation-table.csv'];
    
    % Write the list to file so the IDs can be translated by a bash function
    writetable(table(idList), fullfile(dbDir, 'tmp.csv'), 'WriteVariableNames', false)
    [s, output] = unix([fullfile(dbDir, 'mapIDs.sh'), ' ',...
        fullfile(dbDir, 'tmp.csv'),' ',...
        num2str(targetID), ' ',...
        dbFile]);
    if s
        error('The mapping did not work, most likely the translating file does not exist or is in wrong format')
    end
    
    % Process the terminal output (last line will be always empty
    trList = regexp(output(1:end-1), '\n', 'split');
end

if verbose
    fprintf('\nTranslated %3.2f%%\n\n', 100*sum(~cellfun('isempty', trList))/numel(idList));
end
trList = reshape(strtrim(trList), numel(trList), 1);
trList = regexprep(trList, ';', '|');

end
