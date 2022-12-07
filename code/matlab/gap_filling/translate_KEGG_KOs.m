% get sequence similarity weights from KEGG orthologies
options; clear

fid = fopen('data/RAVEN_dataDir/selfmade/prok90_kegg90/keggdb/ko', 'r');
line = fgetl(fid);
ENTRY = cell(22767, 1);
DBLINKS = cell(22767, 1);
i = 0;
while ischar(line)
    
    ko = contains(line, 'ENTRY');
    
    if any(ko)
        i = i + 1;
        ko = regexp(line, 'K[0-9]{5}', 'match');
        ko = char(ko);
        fprintf('\n%i:\t%s', i, ko)
        ENTRY{i} = ko;
        DBLINKS(i) = {''};
        ko = false;
    end
    
    
    link = contains(line, 'DBLINKS');
    
    if any(link)
        link = regexp(line, 'R[0-9]{5}', 'match');
        link = strjoin(link, '|');
        fprintf('\t%s', link)
        DBLINKS{i} = link;
        link = false;
    end
        
    line = fgetl(fid);
end
fclose(fid);

koTable = cell2table([ENTRY, DBLINKS], 'VariableNames', {'ENTRY', 'RXNS'});

% Translate KEGG IDs to MNXref IDs
for i=1:size(koTable, 1)
    if ~isempty(koTable.RXNS{i})
        rxns = strsplit(koTable.RXNS{i}, '|');
        rxns = cellfun(@(x)translateIDs(x, 'rxn', [], 'KEGG', 'MNXref', false),...
            rxns, 'UniformOutput', false);
        rxns = cellfun(@char, rxns, 'UniformOutput', false);
        rxns = strjoin(rxns, '|');
        koTable.RXNS{i} = strcat('|', rxns, '|');
    end
end

save('data/gap-filling/sequence-similarity/translation-KO-MNXref.mat', 'koTable')
