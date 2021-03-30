%% construct a DAG from the ChEBI ontology
options

% read chebi ontology file
ontologyTab = readtable(ontologyFile, 'ReadVariableNames', true,...
    'Delimiter', '\t');

% empty sparse matrix dimension terms x terms
nz = regexp(ontologyTab.is_a, '\d+');
nz = numel([nz{:}]');
n = size(ontologyTab,1);
ontologyMat = spalloc(n,n,nz);
 
% for every term, set the entries in the columns of its parents to 1
for i=1:n
    
    if mod(i,5000)==0
        disp(i)
    end
    
    row = zeros(1, n);
    
    if ~isempty(ontologyTab.is_a{i})
        
        % get the parent terms for the current term
        parents = regexp(ontologyTab.is_a(i), '\|', 'split');
        parents = [parents{:}]';
        
        % find the row number of its parents in the table
        parent_idx = cellfun(@(x)find(ismember(ontologyTab.ID, str2double(x))), parents,...
            'UniformOutput', false);
        parent_idx = [parent_idx{:}]';
        
        % set the columns at the parent indices to 1 for the current row
        row(parent_idx) = 1;
    end
    ontologyMat(i,:) = sparse(row);
end
terms = ontologyTab.ID;
clear ontologyTab row parents parent_idx nz

save(chebiOntologyWS, 'ontologyMat', 'terms')