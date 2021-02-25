function [scores, genes] = sequenceSimilarity(hmmDir, fastaFile, rxnList, transTab, verbose)
% Performs hmmsearch using hmmer-3.2.1 for all Hidden Markow models (HMMs)
% that are available for a list of reactions. The KEGG orthologies (KO) are linked to
% MNXref reactions so the E-values can be directly assgined to the
% reactions. The scores can be directly used as weights in a LP.
%
% Input:
%       char hmmDir:            path to the folder containing the HMMs
%       char fastaFile:         path to the fasta file containing the amino
%                               acid sequences for the ORFs of the organism
%       cell rxnList:           array with the reactions that should be
%                               checked for
%       table transtab:         table containing the columns 'ENTRY' (KEGG
%                               orthology IDs) and 'RXNS' (MNXref reactions
%                               separated by '|': '|MNXR1|MNXR2|*|')
%       logical verbose:        if false, function will remain silent
%                               (default: true)
% Output:
%       double scores:          E-values from hmmsearch for all reactions
%                               in rxnList
%       cell genes:             genes associated with the respective
%                               reactions (only best hit)

if nargin < 4
    error('USAGE: SequenceSimilarity = (hmmDir, fastaFile, rxnList, transTab)')
elseif nargin < 5
    verbose = true;
elseif ~ischar(hmmDir) || ~isdir(hmmDir)
    error('Please check the given hmmDir argument.')
elseif ~ischar(fastaFile) || ~isfile(fastaFile)
    error('Please check the given fastaFile argument.')
elseif ~iscell(rxnList)
    error('Please check the given rxnList argument.')
elseif ~istable(transTab)
    error('Please check the given transTab argument.')
end

% [test_hmmer, ~] = unix('hmmsearch');
% if test_hmmer == 127
%     error('Please make sure, the hmmsearch function is on your search path.')
% end

n = numel(rxnList);
% initialize the scores vector
scores = ones(n , 1);
% initialize the genes array
genes = repmat({''}, n, 1);

% reaction ids possibly contain compartment identifiers
rxnList = regexp(rxnList, 'MNXR[0-9]*', 'match');

% name of all KEGG orthology (KO) HMM files
koFiles = dir(hmmDir);
KO = regexp({koFiles.name}, 'K[0-9]{5}', 'match');
KO = cellfun(@char, KO, 'UniformOutput', false);

% decrease the size of the reaction list by setting the entries that are
% not contained in the KO translation to empty
ko_rxns = regexp(transTab.RXNS, 'MNXR[0-9]*', 'match');
ko_rxns = [ko_rxns{:}];
idx_no_match = ~cell2mat(cellfun(@(x)ismember(x, ko_rxns),...
    rxnList, 'UniformOutput', false));

% add '|' to both ends ot the reaction to simplify the search
rxnList = cellfun(@char, rxnList, 'UniformOutput', false);
rxnList = strcat('|', rxnList, '|');
rxnList(idx_no_match) = {''};

for i=1:n
    if ~isempty(rxnList{i})
        if verbose
            fprintf('\n%s (%d/%d)\n',...
                rxnList{i}, i, sum(~cellfun('isempty', rxnList)))
        end
        % matching indices
        idx = contains(transTab.RXNS, rxnList{i});
        % KO(s) associated with the current reaction
        current_KOs = transTab.ENTRY(idx);
        % take only these KOs with a HMM
        current_KOs = intersect(current_KOs, KO);
        
        if ~isempty(current_KOs)
            if verbose
                fprintf('\t==> Performing hmmsearch...\n')
            end
            for j=1:numel(current_KOs)
                % execute hmmsearch for the current KO
                cmd = ['hmmsearch ', fullfile(hmmDir, strcat(current_KOs{j}, '.hmm')), ' ',...
                    fastaFile];
                [status, output] = unix(cmd);
                if contains(output, '[No hits detected that satisfy reporting thresholds]')
                    if verbose
                        warning('No hits found')
                    end
                elseif status == 0
                    % subsect the output for the desired part
                    out_extract = strtrim(extractBetween(output, '-#dom-', 'Domain'));
                    out_split = strsplit(char(out_extract), '\n');
                    lineOfInterest = strtrim(out_split{3});
                    split = strsplit(lineOfInterest, ' ');
                    
                    % E-value for the best hit
                    evalue = str2double(split{1});
                    
                    % gene associated with best hit
                    gene = split{end};
                    
                    % in case there are muliple KOs associated with one reaction, take
                    % the best overall hit
                    if evalue < scores(i)
                        scores(i) = evalue;
                        genes{i} = gene;
                        clear evalue gene
                    end
                else
                    if verbose
                        warning('hmmsearch has non-zero exit status for %s', current_KOs{j})
                    end
                end
            end
            if verbose
                fprintf('done.\n')
            end
        else
            if verbose
                fprintf('\t==> No HMM found\n')
            end
        end
    end
end
end