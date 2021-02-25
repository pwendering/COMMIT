function changed_model = addMetFormulae(model, formulaeTab)
% get a formula for each meabolite
% Input:
%           struct model:           metabolic model with fields 'mets'
%           table formulaeTab:      translation from metabolite IDs to
%                                   formulae (first column IDs, second
%                                   formuale)
% Output:
%           struct changed_model:   additional field 'metFormulas'y
dbIDs = formulaeTab.(1);
dbForm = formulaeTab.(2);
clear formulaeTab

% remove the compartment information
mets = strtok(model.mets, '\[');
met_form = repmat({''}, numel(mets), 1);

for i=1:numel(mets)
    idx = strcmp(mets(i), dbIDs);
    res = dbForm(idx);
    if ~isempty(res)
        met_form(i) = res;
    end
end

model.metFormulas = met_form;
changed_model = model;
end