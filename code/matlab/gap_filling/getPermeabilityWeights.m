function p = getPermeabilityWeights(mets)
% Get a binary vector that contains 1 if a metabolite has a high
% probability of crossing a membrane and 0 otherwise
% Input:
%           cell mets:          array of metabolite IDs (MNXref namespace)
% Output:
%           double p:           1 if diffusion is likely, 0 otherwise

options; clearvars -except topDir

% load the table containing the properties
load(fullfile(topDir, 'data/gap-filling/molecular-properties/properties_MNXref_mets.mat'));

% remove the compartment identifier
mets = strtok(mets, '[');
% subsect the table
ID = property_table.ID;
property_table = property_table(ismember(ID, mets), 2:end);
% predict the permeability
p = zeros(numel(mets),1);
permeability = predictPermeability(property_table);
c=0;
uniq_mets = unique(mets, 'stable');
for i=1:numel(uniq_mets)
    if ismember(uniq_mets(i),ID)
        c=c+1;
        p(strcmp(mets(i), mets)) = permeability(c);
    end
end
end