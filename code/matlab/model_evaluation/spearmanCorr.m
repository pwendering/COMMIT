function p = spearmanCorr(v_1, v_2)
% Calculate Spearman's rank correlation for vetors v_1 and v_2

if ~isequal(numel(v_1),numel(v_2))
    error('The given vectors do not have the same length.')
end

m = numel(v_1);

ranks_1 = floor(tiedrank(v_1));
ranks_2 = floor(tiedrank(v_2));
d = abs(ranks_1 - ranks_2);
p = 1 - ( 6 * sum(d.^2) ) / (m*( m^2 - 1 ));

end