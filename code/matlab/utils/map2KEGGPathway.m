function pathways = map2KEGGPathway(metList, pathwayFile)
% First translates MNXref identifiers to KEGG identifiers and maps these to
% KEGG pathways (ko).
% Input:
%       cellstr metList:        array containing the metabolite IDs
%       char pathwayFile:       path to the parsed pathway.gz file obtained
%                               from KEGG 
%                               (pathway ID \t name \t compounds 
%                               (| - separated))
% Output:
%       cell pathways:          pathways mapping to every metabolite

%% ~~~~~ check input ~~~~~ %%
if ~iscellstr(metList)
    error('metList is not of type cellstr')
elseif ~ischar(pathwayFile)
    error('pathwayFile is not of type char')
end

%% ~~~~~ translate IDs to KEGG ~~~~~ %%
% if compartment identifiers are present, they are removed
metList = translateIDs(strtok(metList, '['), 'met', [], 'MNXref', 'KEGG', false);
metList(cellfun(@isempty, metList)) = {'NA'};
%% ~~~~~ find pathways ~~~~~ %%
% read the file as a table
pTable = readtable(pathwayFile, 'ReadVariableNames', true, 'Delimiter', '\t');
% find the matching indices for every metabolite and associated pathway IDs
pathways = cellfun(@(x)pTable.name(contains(pTable.compounds, strsplit(x, '|'))), metList, 'un', 0);

end
