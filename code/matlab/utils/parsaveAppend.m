function parsaveAppend(filename, model)
% Save a variable in a parfor loop
% Input:
%       char filename:          Name of workspace file
%       char model:             Name of variable

save(filename, 'model', '-append')