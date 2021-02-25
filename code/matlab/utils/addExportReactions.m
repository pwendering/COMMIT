function model = addExportReactions(model, metList)
% Add transport reactions from the cytosol to extracellular space
% Input:
%           struct model:               metabolic model
%           cellstr metList:            list of metabolites
% Output:
%           struct model:               model with additional transport
%                                       reactions

comp1 = '[c]';
comp2 = '[e]';
% remove compartment information
metList = strtok(metList, '[');

% add transport reactions
for i=1:numel(metList)
    model = addReaction(model, strcat('export_', metList{i}),...
        'reactionName', strcat(metList{i}, '_transport_by_diffusion'),...
        'metaboliteList', {strcat(metList{i}, comp1),...
            strcat(metList{i}, comp2)},...
        'stoichCoeffList', [-1 1],...
        'reversible', false,...
        'printLevel', 0);
    
end

end