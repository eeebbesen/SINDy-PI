close all;clear all; clc;
exp = 2;
% exp1
if exp==1
exp1res = load("Results\filter_results.mat");
numReps = 10;
numMethods = 3;

RMSE = zeros(numReps, numMethods);  % pre-allocate

for i = 1:numReps
    RMSE(i,1) = sqrt(mean(exp1res.Results(i).sin.^2));
    RMSE(i,2) = sqrt(mean(exp1res.Results(i).reg.^2));
    RMSE(i,3) = sqrt(mean(exp1res.Results(i).noi.^2));
end

% Flatten data for anovan
RMSE_vec = RMSE(:);                     % all observations in a vector
Method = repmat({'SINDy','Regression','Unfiltered'}, numReps, 1);
Method = Method(:);                     % matching factor vector

% Run one-way ANOVA with replication
[p, tbl, stats] = anovan(RMSE_vec, {Method}, 'model','linear', 'varnames', {'Method'}, 'display', 'on');

% Post-hoc pairwise comparisons
c = multcompare(stats, 'CType','bonferroni')

% Compute residuals: observation minus group mean
groupMeans = grpstats(RMSE_vec, Method, 'mean');  % group means
residuals_vec = RMSE_vec - repelem(groupMeans, numReps);  % subtract group mean from each observation

% Plot histogram and Q-Q plot
figure;
subplot(1,2,1);
histogram(residuals_vec);
title('Histogram of residuals');
xlabel('Residual');
ylabel('Frequency');

subplot(1,2,2);
qqplot(residuals_vec);
title('Q-Q plot of residuals');
% 
% % Jarque-Bera test
% [h,p] = jbtest(residuals_vec);
% if h == 0
%     fprintf('Residuals are roughly normal (fail to reject H0), p = %.3f\n', p);
% else
%     fprintf('Residuals are NOT normal (reject H0), p = %.3f\n', p);
% end
end
if exp==2
%% Load data
allRes = load("Results\all_results_no_0.mat"); 
data = allRes.Results;  % struct with fields: freq, input, noise, trial, result

firstVals = arrayfun(@(s) s.result.RMSE(1), data);
secondVals = arrayfun(@(s) s.result.RMSE(1), data);
combinedRMSE = sqrt(firstVals.^2 + secondVals.^2);

formated  = [[data.trial].',[data.freq].',[data.noise].',[data.input].',combinedRMSE.'];

% Columns in 'formated': [trial, freq, noise, input, RMSE]
trials = formated(:,1);       % replicate number
freq   = formated(:,2);       % factor A
noise  = formated(:,3);       % factor B
input  = formated(:,4);       % factor C
RMSE   = formated(:,5);       % dependent variable

% Combine factors into a cell array for anovan
factors = {freq, noise, input};

% Factor names
factorNames = {'Frequency', 'Noise', 'Input'};

%% --- Apply square-root transformation ---
RMSE_reci = log(RMSE);

%% --- Run three-way factorial ANOVA with transformed data ---
[p, tbl, stats] = anovan(RMSE, factors, ...
                         'model', 'full', ...          % include main effects + all interactions
                         'varnames', factorNames, ...
                         'random', 1, ...              % if 'trial' is random
                         'display', 'on');

% Multiple comparisons (example for first factor)
c = multcompare(stats, 'Dimension', 1, 'CType', 'bonferroni');

%% --- Normality check of residuals ---
residuals_vec = stats.resid;  % residuals after transformation

% Plot histogram and Q-Q plot
figure;
subplot(1,2,1);
histogram(residuals_vec);
title('Histogram of residuals (sqrt RMSE)');
xlabel('Residual');
ylabel('Frequency');

subplot(1,2,2);
qqplot(residuals_vec);
title('Q-Q plot of residuals (sqrt RMSE)');

[h,p] = jbtest(residuals_vec);
if h == 0
    fprintf('Residuals are roughly normal (fail to reject H0), p = %.3f\n', p);
else
    fprintf('Residuals are NOT normal (reject H0), p = %.3f\n', p);
end

figure;
boxplot(residuals_vec);
title('Boxplot of residuals (sqrt RMSE)');
ylabel('Residual');
end
% %% add duplicates
% i = 1;
% while i <= length(data)
%     if data(i).noise == 0
%         % Number of replicas to add
%         numReplicas = 4;
% 
%         % Create replicas
%         replicas = repmat(data(i), 1, numReplicas);
%         for r = 1:numReplicas
%             replicas(r).trial = data(i).trial + r;  % increment trial number
%         end
% 
%         % Insert replicas after current element
%         data = [data(1:i), replicas, data(i+1:end)];
% 
%         % Skip past original + replicas
%         i = i + 1 + numReplicas;
%     else
%         i = i + 1;
%     end
% end
% Results = data;
% 
% % Save modified data back to MAT file
% save("Results\all_results_with_replicas.mat", "Results");
% 
% %% remove 0 noise
% % Assuming 'data' is your struct array
% 
% % Find indices where noise is NOT zero
% keepIdx = [data.noise] ~= 0;
% 
% % Keep only those elements
% data = data(keepIdx);
% 
% % Save the cleaned data back to a MAT file
% Results = data;
% save("Results\all_results_no_0.mat", "Results");