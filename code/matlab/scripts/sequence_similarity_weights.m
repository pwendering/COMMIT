% get sequence similarity weights
disp('-------------------------------------------------------------------')
disp('START')
disp('-------------------------------------------------------------------')

% number of genomes/models
n = 432;

% top folders for fasta file
topDir = '/net/calc1/srv/wendering/genomes/ORFs-aminoacid';
% folders = dir(topDir);
% habitats = regexp({folders.name}, '\w+', 'match');
% habitats = [habitats{:}];
habitats = {'Soil', 'Root', 'Leaf'};
% load the gap-filling database
load('/stud/wendering/Masterthesis/DATA/Gap-filling/database/Universal-model-MNXref-balanced.mat')
% load the translation file from KEGG KOs to MNXref reactions
load('/stud/wendering/Masterthesis/DATA/Gap-filling/sequence-similarity/translation-KO-MNXref.mat')

% initialize the matrix where the weights are stored
seq_sim_mat = zeros(numel(dbModel_MNXref_balanced.rxns), n);
col_labels = repmat({''}, n, 1);
genes = cell(n, 1);
% save('/stud/wendering/Masterthesis/DATA/Gap-filling/sequence-similarity/sequence_similarity.mat',...
%     'seq_sim_mat', 'col_labels', 'genes')
load('/stud/wendering/Masterthesis/DATA/Gap-filling/sequence-similarity/sequence_similarity.mat')
% directory where the Hidden Markow models for each KEGG KO are stored
hmmDir = '/stud/wendering/Masterthesis/DATA/Gap-filling/HMMs-KEGG-orthologies';

c = 0;
for i=1:numel(habitats)
    disp(habitats{i})
    disp('------')
    % get fasta file names
    fastaDir = fullfile(topDir, habitats{i});
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
            save('/stud/wendering/Masterthesis/DATA/Gap-filling/sequence-similarity/sequence_similarity.mat',...
                'seq_sim_mat', 'col_labels', 'genes', '-append')
        end
    end
    disp('-------------------------------------------------------------------')
    
end
disp('-------------------------------------------------------------------')

        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        