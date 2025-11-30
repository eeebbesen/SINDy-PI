%% Close all, clear all, clc
close all;clear all; clc;
[status,msg] = mkdir('Results');
addpath('Functions')
set(0,'defaulttextInterpreter','latex')
format long
%%
% Hz, u, u_test, noise
fprintf('Run started at %s\n', datestr(now,'HH:MM:SS'));
m1=0.002000470448860; m2=0.002153227319931;

filterResultsFile = 'Results/filter_results.mat';
Results = [];   % initialize empty

% Pendulum arm length
l1=0.0925; l2=0.063;

% Damping ratio
Jm1=0.000095751475154; Jm2=0.000899016338872;
b1=0.005759843843316; b2=0.012080972309640;
tau1=1e-12; tau2=0.012080972309640;
tanh_k=1000;

Control = 1;
Shuffle = 0;

freq = 1000;
noise  = 0;

dataFiles = { ...
    'C:\Users\ppiuq\Desktop\Test data 9\desync chirp f1.mat'
    % add more here if needed
};

data_range = 4500:20500;
data_type_range = 13:14;
step = 1;

Data = [];
u = [];
tspan = [];

%% Get training data
for f = 1:length(dataFiles)
    loaded = load(dataFiles{f});
    vars = fieldnames(loaded);
    Y = loaded.(vars{1}).Y;
    X = loaded.(vars{1}).X;

    % Extract states, input, time
    stateArray = arrayfun(@(x) x.Data(1,data_range), Y(data_type_range), 'UniformOutput', false);
    uArray = arrayfun(@(x) x.Data(1,data_range), Y(1:2), 'UniformOutput', false);
    tArray = arrayfun(@(x) x.Data(1,data_range), X(1), 'UniformOutput', false);

    Data_tmp = [stateArray{1}.', stateArray{2}.'];
    u_tmp = [uArray{1}; uArray{2}].*-0.33;
    u_tmp(1,:) = u_tmp(1,:)/1.5;
    tspan_tmp = tArray{1};

    % Downsample
    idx = 1:step:length(tspan_tmp);
    Data_tmp = Data_tmp(idx, :);
    u_tmp = u_tmp(:, idx);
    tspan_tmp = tspan_tmp(idx);

    % Adjust time so datasets concatenate smoothly
    if isempty(tspan)
        t_offset = 0;
    else
        t_offset = tspan(end) + (tspan_tmp(2)-tspan_tmp(1));
    end
    tspan_tmp = tspan_tmp - tspan_tmp(1) + t_offset;

    % Concatenate
    Data = [Data; Data_tmp];
    u = [u, u_tmp];
    tspan = [tspan; tspan_tmp];
end

%% get validation data
data_range = 5000:7000;
data_type_range = 13:14;
step = 1;

loaded_test = load('C:\Users\ppiuq\Desktop\Test data 9\constant freq f1.mat');
vars_test = fieldnames(loaded_test);
Y_test = loaded_test.(vars_test{1}).Y;         
X_test = loaded_test.(vars_test{1}).X;

stateArray_test = arrayfun(@(x) x.Data(1,data_range), Y_test(data_type_range), 'UniformOutput', false);
uArray_test = arrayfun(@(x) x.Data(1,data_range), Y_test(1:2), 'UniformOutput', false);
tArray_test = arrayfun(@(x) x.Data(1,data_range), X_test(1), 'UniformOutput', false);

Data_test = [stateArray_test{1}.', stateArray_test{2}.'];
u_test = [uArray_test{1}; uArray_test{2}].*-0.33;
u_test(1,:) = u_test(1,:)/1.5;
tspan_test = tArray_test{1};

idx_test = 1:step:length(tspan_test);
Data_test = Data_test(idx_test, :);
u_test = u_test(:, idx_test);
tspan_test = tspan_test(idx_test);

%% get testing data

data_range = 30150:155150;
data_type_range = [5,6,13,14];
step = 5;

loaded_val = load('C:\Users\ppiuq\Desktop\Test data 2\m03 sin .1t 25f.mat');
loaded_val = load('C:\Users\ppiuq\Desktop\Test data 2\m01 0.1t 25f.mat');
vars_val = fieldnames(loaded_val);
Y_val = loaded_val.(vars_val{1}).Y;         
X_val = loaded_val.(vars_val{1}).X;

stateArray_val = arrayfun(@(x) x.Data(1,data_range), Y_val(data_type_range), 'UniformOutput', false);
uArray_val = arrayfun(@(x) x.Data(1,data_range), Y_val(1:2), 'UniformOutput', false);
tArray_val = arrayfun(@(x) x.Data(1,data_range), X_val(1), 'UniformOutput', false);

Data_val = [stateArray_val{3}.', stateArray_val{4}.'];
Data_val_ref = [stateArray_val{1}.', stateArray_val{2}.'];
u_val = [uArray_val{1}; uArray_val{2}].*-0.33;
u_val(1,:) = u_val(1,:)/1.5;
tspan_val = tArray_val{1};

idx_val = 1:step:length(tspan_val);
Data_val = Data_val(idx_val, :);
Data_val_ref = Data_val_ref(idx_val, :);
u_val = u_val(:, idx_val);
tspan_val = tspan_val(idx_val);

%% Common parameters
dt = tspan(2) - tspan(1);
dData = 0;
dData_test = 0;

    % Define your input signals as a cell array

                %result = run_double_pendulum_plot(u, u_test, 0, Data, dData, Data_test, dData_test, dt, tspan, tspan_test);
                % store in struct
                % Results(counter).freq   = freq;
                % Results(counter).input  = k;
                % Results(counter).noise  = noise;
                % Results(counter).trial  = i;
                % Results(counter).result = result;

                % save after each run (overwrite file with updated Results)
                % 
                % save(resultsFile, 'Results');
                % fprintf('Run %d completed at %s\n', counter, datestr(now,'HH:MM:SS'));
                % 
                % counter = counter + 1;
%%EKF
%28000 length
Data_val1 = Data_val(1:end-1,1);
Data_val2 = Data_val(1:end-1,2);
Data_val_parts1 = reshape(Data_val1, 2500, 10);
Data_val_parts2 = reshape(Data_val2, 2500, 10);

Data_val_ref1 = Data_val_ref(1:end-1,1);
Data_val_ref2 = Data_val_ref(1:end-1,2);
Data_val_ref_parts1 = reshape(Data_val_ref1, 2500, 10);
Data_val_ref_parts2 = reshape(Data_val_ref2, 2500, 10);

u_val1 = u_val(1,1:end-1);
u_val2 = u_val(2,1:end-1);
u_val_parts1 = reshape(u_val1, 2500, 10);
u_val_parts2 = reshape(u_val2, 2500, 10);

tspan_val_trim = tspan_val(1:end-1);
tspan_val_parts = reshape(tspan_val_trim, 2500, 10);

for i=1:size(tspan_val_parts,2)
    Data_itr = [Data_val_parts1(:,i),Data_val_parts2(:,i)];
    tspan_itr = tspan_val_parts(:,i).';
    u_itr = [u_val_parts1(:,i).';u_val_parts2(:,i).'];
    Data_ref_itr = [Data_val_ref_parts1(:,i),Data_val_ref_parts2(:,i)];

    [Data_filtered,meas] = EKF(Data_itr, tspan_itr, u_itr, 0);

    plot(tspan_itr,Data_filtered(1,:))
    legend("sindy")
    errors = Data_ref_itr - Data_filtered(1:2,:).';
    RMSE_sindy = sqrt(mean(errors.^2, 1));
    
    hold on
    [Data_filtered,meas] = EKF(Data_itr, tspan_itr, u_itr, 1);
    plot(tspan_itr,Data_filtered(1,:))
    errors = Data_ref_itr - Data_filtered(1:2,:).';
    RMSE_reg = sqrt(mean(errors.^2, 1));
    
    hold on
    plot(tspan_itr,Data_itr(:,1))
    errors = Data_ref_itr - Data_itr;
    RMSE_noise = sqrt(mean(errors.^2, 1));
    
    hold on
    plot(tspan_itr,Data_ref_itr(:,1))

    Results(i).sin = RMSE_sindy;
    Results(i).reg = RMSE_reg;
    Results(i).noi = RMSE_noise;
end

save(filterResultsFile, 'Results');




