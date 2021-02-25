function [status, fileList] = createModelDir(models, outDir)
% Write the metabolic models from a cell array to single .mat files in a
% given directory
% Input:
%       cell models:        contains metabolic models as struct
%       char dir:           absolute path to directory to write model
%                           files into
% Output:
%       logical status:     whether or not the operation was successful
%       cell fileList:      list of files in the directory

if ~exist('outDir', 'dir')
    status = mkdir(outDir);
else
    status = true;
end



if status
    
    % get model IDs that will be used as file names
    fileList = cellfun(@(x)strtok(x.id, '_'), models,...
        'UniformOutput', false);
    
    for i=1:numel(models)
        model = models{i};
        save(fullfile(outDir, fileList{i}), 'model')
    end
    
else
    
    warning('Directory could not be created')
    fileList = {''};
end

end