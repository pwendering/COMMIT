function norm_tab = CSSNormalization(tab, N)
% Sumulative-sum scaling normalization (CSS)
% Taken from Paulson et al., 2013, Nature Methods

if nargin < 2
    error('USAGE: CSSNormalization(tab, N)')
elseif ~istable(tab) && ~ismatrix(tab)
    error('Argument ''tab'' has to be a table or a matrix')
elseif ~isnumeric(N)
    error('Scaling factor ''N'' has to be numeric')
end

if istable(tab)
    rownames = tab.Properties.RowNames;
    colnames = tab.Properties.VariableNames;
    tab = table2array(tab);
end

% set entries with 1 to 0
tab(tab==1) = 0;

% 95% percentile for each sample (column)
ql = zeros(size(tab, 2), 1);
for i=1:numel(ql)
  ql(i) = quantile(tab(:,i), 0.95);
end


% sum of entries up to ql(j) in each column
sl = zeros(size(tab, 2), 1);
for i=1:numel(ql)
  sl(i) = sum(tab(tab(:,1)<=ql(i), i));
end

% normalization
norm_tab = tab;
for i=1:size(tab,1)
    for j=1:size(tab, 1)
        norm_tab(i,j) = ( tab(i,j) / sl(j) ) * N;
    end
end

norm_tab(isnan(norm_tab)) = 0;

if exist('colnames', 'var')
    norm_tab = array2table(norm_tab,...
        'VariableNames', colnames,...
        'RowNames', rownames);
end

end

