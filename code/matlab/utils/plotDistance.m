function plotDistance(matrix, labels, save, filename, plot_title)
% Plot the distance matrix as a heatmap
% Input:
%           double matrix:          distance matrix
%           cell labels:            labels for the hatmap plot
%           logical save:           if true, file will be saved in the 
%                                   filename path with the given file extension (default false)
%           char filename:          path to save to file
%           char plot_title:        title for the plot

if nargin < 5 || isempty(filename) || ~ischar(filename) || ~islogical(save) || ~ischar(plot_title)
    
    save = false;
    plot_title = '';
    fprintf('\nEither you did not give a filename, it is not of type char or "save" is not of logical type, your plot will not be saved\n')
end

if nargin >= 2
f = figure('visible','off');
h = heatmap(labels, labels, matrix, 'Title', plot_title);
h.ColorMethod = 'none';
colormap(h, flipud(gray(7)))

if save
    saveas(f, filename)
else
    f.Visible = 'on';
end

end
    