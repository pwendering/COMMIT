function [sampleID, taxonomy, keggID] = getKEGGTaxID(wgsTaxonomyFile)
%% [sampleID, taxonomy, keggID] = getKEGGTaxID(wgsTaxonomyFile)
% Extracts the taxonomy of the csv file at the genus level,
% extracts the according organism ids using the KEGG REST API and selects a
% species arbitrarily. If the resolution onyl reaches to family level, the
% output will be empty for this sample.
%
% Input:
%           wgsTaxonomyFile:        csv file that has one leftmost column
%                                   of the sample names and descending
%                                   taxonomy units in the following columns
%                                   must be given as an absolute path
%
% Output:
%           cell sampleID:          cell array of sample ids
%           cell taxonomy:          cell array of the most precise taxonomy
%                                   units
%           cell keggID:            cell array of the selected kegg
%                                   organism ids

% Check if an input was given
if nargin<1
    error('No Taxonomy file given');
end

% Try to open the input file
try 
    fid = fopen(wgsTaxonomyFile);
    line = fgetl(fid);
catch
    error(['Cannot open the TaxonomFile you provided at ', wgsTaxonomyFile]);
end

% Declare the cell arrays to return
sampleID = cell(0,1);
taxonomy = cell(0,1);
keggID = cell(0,1);

% write output to file
name = strsplit(wgsTaxonomyFile, '/');
path = strsplit(wgsTaxonomyFile, name{end});
path = path{1};

fid_out = fopen(strcat(path, 'KEGGTaxonomyIDs.csv'), 'w');

% Read the lines in the taxonomy file; first line is header
line = fgetl(fid);
while ischar(line)
    split = cellstr(strsplit(line, '\t'));
    sampleID{end+1} = split{1};
    taxonomy{end+1} = split{end};
    % use the KEGG API to get the organism ids
    if ~isempty(taxonomy{end})
        [status, output] = unix(['wget -q -O - http://rest.kegg.jp/find/genomes/',...
            taxonomy{end}]);
        if status == 0 && ~isempty(output)
            output = strsplit(strtrim(output), '\n');
            % the organism name should not contain 'sp.'
            output = output(~contains(output, 'sp.'));
            output = output(contains(output, taxonomy{end}));
            output = cellfun(@(x)strsplit(x,{'\t',','}), output, 'UniformOutput', false);
            ids = cell(0,1);
            for i=1:numel(output)
                if numel(output{i})==3
                    ids(end+1) = output{i}(2);
                end
            end
            % sample one id at random
            if ~isempty(ids)
                keggID(end+1) = datasample(ids, 1);
            end
        else
            keggID(end+1) = {''};
        end
    else
        keggID(end+1) = {''};
    end
    fprintf(fid_out,'%s,%s,%s\n', sampleID{end}, taxonomy{end}, keggID{end});
    line = fgetl(fid);
end
fclose(fid);
fclose(fid_out);
end

