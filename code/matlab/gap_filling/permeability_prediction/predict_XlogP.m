%% Read Input
options; clear
% table with decriptors
property_table = readtable('data/gap-filling/molecular-properties/METS-PROPERTIES.csv',...
    'ReadVariableNames', true);
property_table = property_table(:,[1 5 2:4 6:end]);

% select rows that have a XlogP value different from 'NULL'
xlogP_set = property_table(cellfun(@(x)~isequal(x, 'NULL'), property_table.xlogp), :);
% cast character to number array
xlogP_set.xlogp = str2num(char(xlogP_set.xlogp));


%% Split data into training, testing and validation set

% set aside validation set
n_validation = floor(0.05*size(xlogP_set,1));
idx_validation = randsample(size(xlogP_set,1), n_validation);
validation_set = xlogP_set(idx_validation, :);

% exclude the validation set
S = xlogP_set;
S(idx_validation,:) = [];


% Training set
n_train = floor(0.8*size(S,1));
training_set = xlogP_set(randsample(size(S,1), n_train), :);

% Testing set (excluding training set IDs)
testing_set = xlogP_set(~ismember(S.ID, training_set.ID), :);

% create the table for the training set last column is the response
training_set = training_set(:,[3:end,2]);
testing_set = testing_set(:, [3:end,2]);

%% Investigation of the training set
labels = training_set.Properties.VariableNames;
M = table2array(training_set);
% Correlation between the descriptors
correlationPlot(M, labels);
saveas(gcf, 'figures/xlogP_prediction/correlation_features.jpg');

% Histograms showing the desitibutions of each descriptor together with
% XlogP
FaceColor = [.1 0.1 .6];
m=3; n=3; 
pos = 1:numel(labels);
for i=1:numel(labels)
    % position
    subplot(m,n,pos(i))
    % histogram plot
    histogram(M(:,i))
    % fit of normal distibution to data
    h = histfit(M(:,i));
    % select the axis to obtain x and y limits
    ax = gca;
    % label
    text(0.4*ax.XLim(2), 0.7*ax.YLim(2), labels{i}, 'FontWeight', 'bold', 'FontSize', 18)
    % p-value of Kolmogorov-Smirnov test
%     [~,pval] = kstest(M(:,i));
%     text(0.8*ax.XLim(2), 0.9*ax.YLim(2), strcat('p=',num2str(pval, 3)))
    % set the values for appearance
    set(h(1), 'facecolor', FaceColor)
    set(h(2), 'color', 'r', 'LineStyle', ':')
end
% suptitle('Histograms showing the distribution of the descriptors')
% save figure
saveas(gcf, 'figures/xlogP_prediction/histograms_features.jpg');




%% testing set
% Correlation showed that hbondacc, complexity and heavycnt can be excluded
subset_idx = [1 2 5 7 8];
training_set = training_set(:, subset_idx);
testing_set = testing_set(:, subset_idx);



%% Linear regression
% complexity has no significant influence on the output of the model
% (ANOVA)
nbins = 100;
FaceColor = [.1 0.1 .6];
% Fit the linear model
lm = fitlm(training_set);
% predictors and response for the testing set
X_test = testing_set(:,1:end-1);
Y_true = testing_set.xlogp;
% predict XlogP for the test set properties
Y_pred = predict(lm, X_test);
% calculate squared errors
sqe_lm = ( Y_true - Y_pred ).^2;
% Calculate the adjusted R-squared
SS_tot = sum(( Y_true - mean(Y_true) ).^2);
SS_res = sum(( Y_true - Y_pred ).^2);
rsq = 1 - SS_res / SS_tot;
n_sample = size(X_test,1);
p = size(X_test,2);
Rsq_adj = 1-(1-rsq)*( n_sample - 1) / ( n_sample - p - 1 );

% Plot squared errors on histogram and save the figure
histogram(sqe_lm, nbins,...
    'FaceColor', FaceColor,...
    'FaceAlpha', 0.7...
    )
legend({strcat('linear model (max=', num2str(max(sqe_lm), 2),...
    '; rmse=', num2str(sqrt(mean(sqe_lm)), 2), ')')}, 'FontSize', 12)
xlabel('squared error')
ylabel('count')
ax = gca;
text(.5*ax.XLim(2), .5*ax.YLim(2), strcat('R^2_{adj}=',...
    num2str(Rsq_adj,2)), 'FontSize', 16)
% title('Squared errors of the linear model')
ax.TickLength = [.005 .2];
saveas(gcf, 'figures/xlogP_prediction/squared_errors_lm.jpg');



%% K nearest neighbor classifier
N = [1 2 3 4 5 6 7 8 9 10 13 16 20];
% C = {[0 0 0], [0 0 1], [1 0 0], [0 1 0.8], [0.5 0 0.5]};
rmse = [];
SS_tot = zeros(numel(N), 1);
SS_res = zeros(numel(N), 1);
dist_method = 'seuclidean'; % euclidean distance scaled by sd

% select predictor values of the training set
X_train = table2array(training_set(:,1:end-1));
% response of training set
Y_train = training_set.xlogp;
% predictors of the testing set
X = table2array(testing_set(:,1:end-1));
% true XlogP value
Y_true = testing_set.xlogp;
     
for i=1:numel(N)
    n_neighbors = N(i);
    
    % find indices of nearest neighbors
    idx_nn = knnsearch(X_train,...
        X, 'K', n_neighbors, 'Distance', dist_method);
    
    Y = mean(Y_train(idx_nn), 2);
    SS_tot(i) = sum(( Y_true - mean(Y_true) ).^2);
    SS_res(i) = sum(( Y_true - Y).^2);
    
    % for sqerr distribution and rmse
    sqe_knn = ( Y_true - Y ).^2;
    rmse(end+1) = sqrt(mean(sqe_knn));
%     histogram(sqe_knn, nbins,...
%         'DisplayStyle', 'stairs',...
%         'EdgeColor', 'black',...
%         'BinLimits', [0 100],...
%         'LineWidth', 0.8 ...
%         )
%     hold on
%     txt{i}= strcat('kNN (N=', num2str(n_neighbors), ';', dist_method, '; max=',...
%         num2str(max(sqe_knn)), '; rmse= ', num2str(rmse(end),2), ', )');
end

% Settings for the histogram plot
% legend(txt)
xlim([0 50])
xlabel('squared error')
ylabel('count')
text(25, 200, strcat('n = ', num2str(size(testing_set,1))), 'FontSize', 16)
title('Histogram of squared errors for the linear and knn regression model')

% Plot RMSE and adjusted R-squared against the number of neighbors
% RMSE
yyaxis left
col = [0.6 0.1 0.2];
plot(N,rmse, 'Color', col)
xlabel('number of neighbors used')
ylabel('RMSE')
ax = gca;
% ax.YLim(1)=0;
ax.YColor = col;
% title('Root mean square error (RMSE) with different N')
hold on

% coefficient of determination (Rsquared)
yyaxis right
rsq = 1 - SS_res./SS_tot;
% adjust Rsq
Rsq_adj = 1-(1-rsq)*( size(X,1) - 1) / ( size(X,1) - size(X,2) - 1 );
col = 'b';
plot(N,Rsq_adj, 'Color', col)
xlabel('number of neighbors used')
ylabel('R^2_{adj}', 'Color', 'black')
ax = gca;
ax.YColor = col;
ax.YLim(1)=0.8;
ax.TickLength = [.005 .2];
% title('RMSE and adjusted R-squared values with different N')
rectangle('Position', [0.5 0 1 ax.YLim(2)], 'FaceColor', [.7 .7 .7 0.3],...
    'LineStyle', 'none', 'Curvature', 0.2);

saveas(gcf,  'figures/xlogP_prediction/RMSE_Rsq_kNN.jpg');

%% Compare linear regression with kNN (N=1)

% kNN 
% knn_fit = fitcknn(X_train, Y_train, 'NumNeighbors', 1,...
%         'Distance', dist_method, 'PredictorNames', labels(2:end),...
%         'ResponseName', labels{2});
idx_nn = knnsearch(X_train,...
        X, 'K', 1, 'Distance', dist_method);
Y = mean(Y_train(idx_nn), 2);
sqe_knn = ( Y_true - Y ).^2;


nbins = 100;
binmax = max([max(sqe_lm), max(sqe_knn)]);
% linear model
histogram(sqe_lm, nbins,...
    'FaceColor', 'b',...
    'FaceAlpha', 0.7,...
    'BinLimits', [0 binmax],...
    'LineWidth', 0.8 ...
    )

hold on

% kNN
histogram(sqe_knn, nbins,...
        'DisplayStyle', 'bar',...
        'FaceColor', [0.6 0.1 0.2],...
        'EdgeColor', 'black',...
        'BinLimits', [0 binmax],...
        'LineWidth', 0.8 ...
        )

legend({'linear model', 'kNN model'}, 'FontSize', 12)
xlabel('Squared error')
ylabel('count')
title('Squared errors to the true response in the testing set',...
    'FontSize', 18)

saveas(gcf,  'figures/xlogP_prediction/squared-errors-lm-knn.jpg');


%% 10-fold cross-validation
n_neighbors = 1;
dist_method = 'seuclidean';

k = 10;
n = size(S,1);

% group size
interval = ceil(n/k);
% RMSE
rmse = zeros(k,1);
% adjusted R_squared
rsq_adj = zeros(k,1);
% shuffle the dataset
S = S(randsample(n,n),:);
disp('-----------------------------------------')
fprintf('%d-fold cross-validation...\n', k)
disp('-----------------------------------------')
for i=1:k
    
    % select current group for testing
    start_idx = (i-1)*interval+1;
    if i<k
        end_idx = i*interval;
    else 
        end_idx = size(S,1);
    end
    testing_set = S(start_idx:end_idx, :);
    
    % the remaining rows are training data
    training_set = S(~ismember(S.ID, testing_set.ID), :);
    
    % create the table for the training set last column is the response
    training_set = training_set(:,[3,4,7,9,2]);
    testing_set = testing_set(:, [3,4,7,9,2]);

    % select predictor values of the training set
    X_train = table2array(training_set(:,1:end-1));
    % response of training set
    Y_train = training_set.xlogp;
    
    % predictors of the testing set
    X_test = table2array(testing_set(:,1:end-1));
    % true XlogP value
    Y_true = testing_set.xlogp;
    
    % make the prediction
    idx_nn = knnsearch(X_train,...
        X_test, 'K', 1, 'Distance', dist_method);
    Y_pred = mean(Y_train(idx_nn), 2);
    sqe_knn = ( Y_pred - Y_true ).^2;
    
    % RMSE
    rmse(i) = sqrt(mean(sqe_knn));
    
    % calculate the adjusted R squared
    SS_tot = sum(( Y_true - mean(Y_true) ).^2);
    SS_res = sum(( Y_true - Y_pred ).^2);
    
    rsq = 1 - ( SS_res / SS_tot );
    n_sample = size(X_test,1); % number of rows in the test set
    p = size(X_test,2); % number predictors
    rsq_adj(i) = 1 - ( 1 - rsq ) * (n_sample - 1) / (n_sample - p - 1 );
    
end
fprintf('Mean R-squared_adj:\t%.2f\n', mean(rsq_adj))
fprintf('Mean RMSE:\t\t%.2f\n', mean(rmse))
disp('-----------------------------------------')

%% Assess quality using validation set
disp('-----------------------------------------')
fprintf('Skill assessment using validation set with the remaining set as training set')
disp('-----------------------------------------')
X_valid = table2array(validation_set(:, [3,4,7,9]));
Y_valid = validation_set.xlogp;

% now use the whole training set
X_train = table2array(S(:, [3,4,7,9]));
Y_train = S.xlogp;

idx_nn = knnsearch(X_train, X_valid, 'K', n_neighbors, 'Distance', dist_method);
Y_pred = mean(Y_train(idx_nn), 2);
sqe_knn = ( Y_pred - Y_valid ).^2;

% RMSE
rmse = sqrt(mean(sqe_knn));

% calculate the adjusted R squared
SS_tot = sum(( Y_valid - mean(Y_valid) ).^2);
SS_res = sum(( Y_valid - Y_pred ).^2);

rsq = 1 - ( SS_res / SS_tot );
n_sample = size(X_valid,1); % number of rows in the test set
p = size(X_valid,2); % number predictors
rsq_adj = 1 - ( 1 - rsq ) * (n_sample - 1) / (n_sample - p - 1 );

fprintf('R-squared_adj:\t%.2f\n', rsq_adj)
fprintf('RMSE:\t\t%.2f\n', rmse)
disp('-----------------------------------------')
%% Prediction of unknown XlogP values

% training set
X_train = table2array(xlogP_set(:, [3,4,7,9]));
Y_train = xlogP_set.xlogp;

% properties for prediction
idx_no_logP = find(cellfun(@(x)isequal(x, 'NULL'), property_table.xlogp));
X_prop = property_table(idx_no_logP, :);
X_prop = table2array(X_prop(:, [3,4,7,9]));
% predict XlogP
idx_nn = knnsearch(X_train, X_prop, 'K', n_neighbors, 'Distance', dist_method);
Y_pred = mean(Y_train(idx_nn), 2);
% fill in missing values

property_table.xlogp(idx_no_logP) = {'NaN'};
property_table.xlogp = str2num(char(property_table.xlogp));
property_table.xlogp(idx_no_logP) = Y_pred;

