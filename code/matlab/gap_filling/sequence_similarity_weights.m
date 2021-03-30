% get sequence similarity weights
options

disp('-------------------------------------------------------------------')
disp('START')
disp('-------------------------------------------------------------------')

% number of genomes/models
n = 432;

orfDir = fullfile(topDir, 'data', 'genomes', 'ORFs-aminoacid');
names = dir(orfDir);
names = regexp({names.name}, '\w+', 'match');
habitats = [names{~cellfun(@isempty,names)}];

% load the gap-filling database
load(dbFile)
% load the translation file from KEGG KOs to MNXref reactions
load(fullfile(topDir, 'data/gap-filling/sequence-similarity/translation-KO-MNXref.mat'))

% initialize the matrix where the weights are stored
seq_sim_mat = zeros(numel(dbModel_MNXref_balanced.rxns), n);
col_labels = repmat({''}, n, 1);
genes = cell(n, 1);

% save workspace already (can be loaded with intermediate results and
% completed later if the script was abortet for some reason)
save(seq_sim_workspace, 'seq_sim_mat', 'col_labels', 'genes')
% load(seq_sim_workspace)

% directory where the Hidden Markow models for each KEGG KO are stored
hmmDir = fullfile(topDir, 'data/gap-filling/HMMs-KEGG-orthologies');

c = 0;
for i=1:numel(habitats)
    disp(habitats{i})
    disp('------')
    
    % get fasta file names
    fastaDir = fullfile(orfDir, habitats{i});
    fastaFiles = dir(fastaDir);
    fastaFiles = regexp({fastaFiles.name}, '\w+\.faa', 'match');
    fastaFiles = [fastaFiles{:}];
    for j=1:numel(fastaFiles)
        c = c + 1;
        fastaFile = fastaFiles{j};
        if ~ismember({strtok(fastaFile, '.')}, col_labels)
            col_labels(c) = {strtok(fastaFile, '.')};
            fprintf('HMM search for %s... ', fastaFile)
            fastaFile = fullfile(fastaDir, fastaFile);
            
            % calculate weights
            tic;
            [seq_sim_mat(:, c), genes{c}] = sequenceSimilarity(hmmDir, fastaFile, dbModel_MNXref_balanced.rxns, koTable, false);
            fprintf('(%.2fm)\n', toc/60)
            disp('------')
            save(seq_sim_workspace, 'seq_sim_mat', 'col_labels', 'genes', '-append')
        end
    end
    disp('-------------------------------------------------------------------')
    
end
disp('-------------------------------------------------------------------')

clear topDir
        