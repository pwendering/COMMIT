function model = getMergedModelKEGGMetaCyc(sampleID, keggID, fastaFile, dataDir, outDir, getModelFromID)
% Generate a merged model from KEGG and MetaCyc drafts. The KEGG
% reconstruction is based on pre-trained Hidden Markov Models (HMMs) that
% were generated based on KEGG v87 or v90. The functions of the RAVEN2.0
% toolbox were used to generate the HMMs.
% Input:
%       char sampleID:          descriptive name for the model
%       char keggID:            If given, an additional model will be
%                               reconstructed that is only based on the
%                               keggID. If empty, this step will be skipped
%       char fastaFile:         absolute path to the FASTA file on which
%                               the reconstruction will be based on
%       char dataDir:           absolute path to the directory that
%                               contains the KEGG FTP files
%                               (genes.pep, compound, ko, compound.inchi,
%                               reaction, reaction.lst,
%                               reaction_mapformula.lst),HMMs, FASTA files
%                               and mulitple alignments
%                               (also see RAVEN function getKEGGModelForOrganism)
%       char outDir:            absolute path to the directory where the
%                               .out files are stored, if empty the files
%                               will stored in a temporary directory
% Output:
%       struct model:           merged model

if ischar(keggID)
    fprintf(['\nA KEGG  organism ID was given, an additional model',...
        'will be reconstructed for %s\n'],keggID)
    getModelFromID = 1;
end

% Get the individual reaconstructions

% HMMer cutOff
t_hmm = 10^-50;

% Generate a model using the proteome as input
fprintf('\nObtaining the draft reconstruction based on the given FASTA file\n\n')
KEGG_model = getKEGGModelForOrganism(...
    sampleID,...
    fastaFile,...
    dataDir,...
    outDir,...
    false, false, false, false,...
    t_hmm,0.8,0.3,-1);

if getModelFromID
    fprintf(['\nObtaining the draft reconstruction according to the given KEGG',...
        'organism identifier\n'])
    % Generate a model only with KEGG organism id
    model_ID = getKEGGModelForOrganism(...
        keggID,...
        [],...
        dataDir,...
        outDir,...
        false, false, false, false,...
        t_hmm,0.8,0.3,-1);
    % Merge the KEGG models
    fprintf('\nMerging the two KEGG models\n')
    KEGG_model = mergeModels({model_ID, KEGG_model});
end

% Generate a model using the MetaCyc database (diamond aligner)
fprintf(['\nGenerating a draft reconstruction from MetaCyc',...
    ' using the diamond aligner\n\n'])
MetaCyc_model = getMetaCycModelForOrganism(sampleID,...
    fastaFile);


% Combine the merged KEGG model with the MetaCyc draft model
fprintf('\nCombining the KEGG and MetaCyc models\n\n')
model = combineMetaCycKEGGModels(MetaCyc_model, KEGG_model);
model.id = strcat(sampleID, '_merged');

end





















