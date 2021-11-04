function p = getPermeabilityWeights(mets)
% Get a binary vector that contains 1 if a metabolite has a high
% probability of crossing a membrane and 0 otherwise
% Input:
%           cell mets:          array of metabolite IDs (MNXref namespace)
% Output:
%           double p:           1 if diffusion is likely, 0 otherwise

% load the table containing the properties
load(fullfile('data','gap-filling','molecular-properties','properties_MNXref_mets.mat'))

% remove the compartment identifier
mets = strtok(mets, '[');
% subsect the table
property_table = property_table(ismember(property_table.ID, mets),:);
uniq_mets = unique(mets, 'stable');

% predict the permeability
p = zeros(numel(mets),1);
permeability = predictPermeability(property_table(:,2:end));

for i=1:numel(uniq_mets)
    met_idx = ismember(property_table.ID,uniq_mets(i));
    if sum(met_idx)>0
        p(strcmp(uniq_mets(i), mets)) = permeability(met_idx);
    end
end
end