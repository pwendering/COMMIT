function models = getMultipleModels_par(genomePath, habitat, dataDir, outDir, taxPath, gapFillDir, getModelFromID, writeToFile, modelDir)

% Iterates through all files for the habitat (subfolder of the directory where the
% genomes are located), generates and merges KEGG and MetaCyc-based
% reconstructions.
% Input:
%           char genomePath:            path to the top directory
%                                       where the genomes are located
%           cell habitat:              array containing the names of the
%                                       subfolders in the genomePath
%           char dataDir:               path to the top directory
%                                       where the trained HMMs are stored
%           char outDir:                if empty, .out files will be be
%                                       generated in a temporary directory 
%                                       and deleted after model generation,
%                                       otherwise the .out files will be 
%                                       stored in the given directory
%           char taxPath:               path to the directory where the 
%                                       files contaning the KEGG organism 
%                                       ids is located of where it
%           char gapFillDir:            path to the directory where the
%                                       models are stored, whoch should be 
%                                       used for gap filling 
%                                       (see fillGaps function)
%           logical getModelFromID:     if true, the KEGG organism ids will
%                                       be used, default false
%           logical writeToFile:        if 1, the merged model will be 
%                                       saved as a .xml file, if 2, the
%                                       merged models will be saved in a
%                                       combined .mat file
%           char modelDir:              specifies the location where the 
%                                       model files should be saved

if ~exist('taxPath', 'var')
    taxPath = 0;
else
    wgsTaxonomyFile = fullfile(taxPath, 'wgs_taxonomy.txt');
end

if ~ischar(gapFillDir)
    gapFillDir = [];
end

if ~ischar(modelDir)
    fprintf("\nNo directory for output models given\n")
    writeToFile = 0;
end

if ~exist('writeToFile', 'var') || ~ismember(writeToFile, 0:2)
    writeToFile = 0;
    modelDir = [];
end


% Retrive a list of all relevant FASTA files
allFastaFiles = dir(fullfile(genomePath, habitat));

% Do not include '.' and '..'
allFastaFiles = allFastaFiles(3:end);
allFastaFiles = {allFastaFiles.name};

% Declare cell array that will be returned
models = cell(numel(allFastaFiles),1);

cores = 8;




parfor i=1:numel(allFastaFiles)
    % The .out files have to be written into separate directories:
    out = strcat(outDir, num2str(i));
    if exist(out, 'dir')
        rmdir(out, 's');
    end
    mkdir(out);
   
    % Name of FASTA file and sample id
    fastaFile = strtrim(allFastaFiles{i});
    sampleID = strsplit(fastaFile, '.');
    sampleID = sampleID{1};
    fastaFile = strcat(sampleID, '.faa');
    
    % Get the merged model from KEGG and MetaCyc databases
    fprintf('\n####################################\n\t%s (sample #%d)\n####################################\n',...
        sampleID, i)
    model = getMergedModelKEGGMetaCyc(sampleID, [],...
        fullfile(genomePath, habitat, fastaFile), dataDir, out, getModelFromID);

    % subSystems field is mostly not only of type cell which can lead to
    % problems with writing and reading
    for field=1:numel(model.subSystems)
        if ischar(model.subSystems{field})
            model.subSystems{field} = {model.subSystems{field}};
        end
    end
    
    % Write the model to a .xml file or add it to a workspace
    try
        if writeToFile ~= 0
            if writeToFile == 1
                fprintf('\nWriting model to .xml file\n')
                exportModel(model, fullfile(modelDir, sampleID), true);
            elseif writeToFile == 2
                fprintf('\nSaving model to workspace\n')
                % Save the model as a variable named as the sampleID
                eval(strcat(sampleID, '=', 'model;'));
                % Append the model to the workspace or create .mat file
                ws = fullfile(modelDir, strcat(habitat,'.mat'));
                if ~exist(workspace, 'file')
                    eval(strcat('parsave("',ws,'", ',sampleID,');'));
                else
                    eval(strcat('parsaveAppend("',ws,',"',sampleID,');'));
                end
            end
        end
    catch
        fprintf('\nFile could not be saved\n')
    end
    % Add the current model to the cell array
    models{i} = model;
    % Delete the .out directories
    rmdir(out,'s');
end



end

    
