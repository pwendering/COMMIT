function trList = translateIDs(idList, id_type, dbTable, source, target, verbose)
%% trList = translateIDs(idList, id_type, dbTable, source, target, verbose)
% Translate metabolite IDs from one namespace to another (ModelSEED, KEGG,
% MetaCyc, BiGG, Rhea, EC, ChEBI).
% Input:
%       cellstr idList:             array containing the ids to be translated
%       char id_type:               either 'rxn' or 'met'
%       table dbTable:              (optional) table that contains the 
%                                   translation of identifiers for the given 
%                                   type (if empty attempts to load table
%                                   from 'COMMIT/data/tables/MNXref')
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

if isempty(dbTable)
    dbDir = 'data/tables/MNXref';
    if isequal(id_type, 'met')
        dbTable = readtable(fullfile(dbDir, 'MNXref-met-translation-table.csv'));
    else
        dbTable = readtable(fullfile(dbDir, 'MNXref-rxn-translation-table.csv'));
    end
end

if ~iscellstr(idList)
    if ischar(idList)
        idList = {idList};
    elseif iscell(idList)
        idList = cellfun(@(x)char(x), idList, 'UniformOutput', false);
    else
        error('The input idList is either not given or not of type cellstr')
    end
end

if ~ischar(source) || ~ischar(target) || any(~ismember(id_type, ['rxn' 'met']))
    error('The source and/or target namespace or type definition is incorrect')
end

% Available namespaces
rxnSources = {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED', 'Rhea', 'EC'};
metSources = {'MNXref', 'KEGG', 'BiGG', 'MetaCyc', 'ModelSEED', 'ChEBI', 'NAMES'};

%% Find matched using the MNXref database

% get the respective column that contains the desired names spaces
if isequal(id_type, 'met')
    sourceID = find(contains(metSources, source));
    targetID = find(contains(metSources, target));
else
    sourceID = find(contains(rxnSources, source));
    targetID = find(contains(rxnSources, target));
end

if isempty(sourceID) || isempty(targetID)
    error('The requested namespace is not available')
end

% initialize translated IDs list
trList = repmat({''}, numel(idList), 1);

% Filter source  and target IDs by empty keys and keys that are definitely not
% contained in the given list
sourceIDs = dbTable.(source);
targetIDs = dbTable.(target);
sourceIDs = sourceIDs(contains(dbTable.(source), idList));
targetIDs = targetIDs(contains(dbTable.(source), idList));
clear dbTable;

% pad source IDs and queries with '|'
sourceIDs = strcat('|', sourceIDs);
sourceIDs = strcat(sourceIDs, '|');
idList = strcat(idList, '|');
idList = strcat('|', idList);

% match IDs
match_idx = cellfun(@(x)find(contains(sourceIDs, x)), idList, 'UniformOutput', false);
empty_idx = cellfun('isempty', match_idx);
trList(~empty_idx) = targetIDs([match_idx{:}]);

if verbose
    fprintf('\nTranslated %3.2f%%\n\n', 100*sum(~empty_idx)/numel(idList));
end

% make it a column vector
trList = reshape(strtrim(trList), numel(trList), 1);

trList = strrep(trList, ';', '|');

end
