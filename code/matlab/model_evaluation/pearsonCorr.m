function r = pearsonCorr(v1, v2)
% Calculate the Pearson correlation coefficient between two vectors v1 and
% v2

if iscell(v1)
    v1 = [v1{:}];
end

if iscell(v2)
    v2 = [v2{:}];
end

if ~(numel(v1)==numel(v2))
    error('Vectors have different lengths')
end

m = numel(v1);

av_1 = 1/m * sum(v1);
av_2 = 1/m * sum(v2);
sd_1 = sqrt( 1/m * sum( (v1 - av_1).^2 ));
sd_2 = sqrt( 1/m * sum( (v2 - av_2).^2 ));
cov = 1/m * ( rowvector(v1 - av_1) * columnVector(v2 - av_2) );
r = cov / ( sd_1 * sd_2 );

end

