function brite = map2KEGGBrite(metList, briteFile)
% First translates MNXref identifiers to KEGG identifiers and maps these to
% KEGG br08001 (with manual extension using KEGG, PubChem and ChEBI).
% Input:
%       cellstr metList:        array containing the metabolite IDs
%       char pathwayFile:       path to the parsed br08001.keg file obtained
%                               from KEGG 
%                               (compound \t name \t descending hierarchy  
%                               (| - separated))
% Output:
%       cell pathways:          highest brite level for every metabolite

%% ~~~~~ check input ~~~~~ %%
if ~iscellstr(metList)
    error('metList is not of type cellstr')
elseif ~ischar(briteFile)
    error('briteFile is not of type char')
elseif isempty(metList)
    warning('Empty metabolite list')
    brite = {{''}};
    return
end

%% ~~~~~ translate IDs to KEGG ~~~~~ %%
% if compartment identifiers are present, they are removed
metList_KEGG = translateIDs(strtok(metList, '['), 'met', [], 'MNXref', 'KEGG', false);
% untranslated IDs are re-filled with MNXref IDs (contained in the
% manual extension of brite)
idx_empty = cellfun(@isempty, metList_KEGG);
metList_KEGG(idx_empty) = metList(idx_empty); clear metList
%% ~~~~~ find brite ~~~~~ %%
% read the file as a table
bTable = readtable(briteFile, 'ReadVariableNames', false, 'Delimiter', '\t');
% find the matching indices for every metabolite and associated pathway IDs
brite = cellfun(@(x)bTable.(3)(contains(bTable.(1), strsplit(x, '|'))), metList_KEGG, 'un', 0);
brite = cellfun(@(x)strtok(x, ';'), brite, 'un', 0);
end