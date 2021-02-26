function F = correlationPlot(matrix, labels)
if nargin < 2
    labels = repmat({''}, size(matrix,2), 1);
end
% Scatterplots of each vector pair
n = size(matrix,2);
text_col = [0.1 0.1 0.1];
text_font = 'Courier';
diagonal_fsize = 14;
bg_col = 'white';
prec_corr = 2;
%% positions of plots on grid
pos_mat = zeros(n);
c=0;
for i=1:n
    for j=1:n
        c = c + 1;
        pos_mat(i,j) = c;
    end
end

%% generate the plot
F = figure;
set(F, 'Color', bg_col)
for i=1:n
    for j=i:n
        
        if i~=j
            
            % Scatter plot
            subplot(n,n,pos_mat(i,j))
            
            scatter(matrix(:,j), matrix(:,i), '.', 'MarkerEdgeColor', 'black')
            set(subplot(n,n,pos_mat(i,j)), 'Color', [0.9 0.9 0.9])
            axis on
            
            % pearson correlation on lower triangle field
            subplot(n,n,pos_mat(j,i))
            t = num2str(pearsonCorr(matrix(:,j), matrix(:,i)), prec_corr);
            adj = 0.1*length(labels(i));
            text(0.5-adj, 0.5, t,...
                'FontWeight', 'bold',...
                'FontName', text_font,...
                'FontSize', 16,...
                'Color', text_col)
            axis off
            
        else
            % print labels on diagonal
            subplot(n,n,pos_mat(i,j))
            adj = 0.2*length(labels(i));
            text(0.5-adj, 0.5, labels(i),...
                'FontName', text_font,...
                'FontSize', diagonal_fsize,...
                'FontWeight', 'bold',...
                'FontSize', 20,...
                'Color', text_col)
            axis off
        end
    end
end
end

