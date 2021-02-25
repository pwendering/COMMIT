function changed_model = translateGeneIDs(model, mapping_file)
% Maps the gene IDs from the model to the gene IDs in
% the second column in the mapping file
% Input:
%       struct model:           model containing the field 'genes'
%       char mapping_file:      path to the mapping file
% Output:
%       struct changed_model:   model with updated gene IDs

% read the mapping file
fid = fopen(mapping_file);
F = textscan(fid, '%s\t%s');
fclose(fid);
OLD = F{1};
NEW = F{2};
clear F fid

% translate the gene IDs
genes = model.genes;
for i=1:numel(genes)
    idx = strcmp(genes(i), OLD);
    res = unique(NEW(idx));
    if numel(res) == 1
        genes(i) = res;
    elseif isempty(res)
        warning('%s has no hits in the amino acid ORF file', genes{i})
        genes(i) = {''};
    else
        warning('%s has multiple hits in the amino acid ORF file', genes{i})
        genes(i) = {''};
    end
end
model.genes = genes;
changed_model = model;
end