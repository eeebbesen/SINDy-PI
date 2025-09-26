@ -0,0 +1,361 @@
function results = run_double_pendulum(Hz, u, u_test, noise, Data, dData, Data_test, dData_test,dt,tspan,tspan_test,state0_test)
% RUN_DOUBLE_PENDULUM  Full workflow for SINDy-PI double pendulum inference
%
% Inputs:
%   Hz       - Sampling frequency (Hz)
%   u        - Control input for training (2 x N)
%   u_test   - Control input for testing (2 x N_test)
%   noise    - Noise level (scalar)
%
% Output:
%   results  - Struct with all generated variables (Data, dData, models, etc.)

    [~,~] = mkdir('Results');
    addpath('Functions')


    % Pendulum mass
    m1=0.002000470448860; m2=0.002153227319931;

    % Pendulum arm length
    l1=0.0925; l2=0.063;

    % Gravity constant and the first pendulum arm length
    g=0; L=0.2667;

    % Damping ratio
    Jm1=0.000095751475154; Jm2=0.000899016338872;
    b1=0.005759843843316; b2=0.012080972309640;
    tau1=1e-12; tau2=0.012080972309640;
    tanh_k=1000;
    params_est = 0;

    %% Define trajectory

    %% Declare settings
    Control=1;
    Partial = 1;
    Smoothing = 0;
    sgolay_order = 4;
    sgolay_window = 45;
    Sindy = 1;
    Integral = 1;
    Shuffle = 0;
    ProvidedData = 0;
    Nonlinlsq = 0;
    Compare = 0;
    saveRegODE = 0;
    saveSindyODE = 0;

    %% Add noise
    

    %% === Smoothing and velocity/acceleration estimation ===
if Partial==1 && Smoothing==0
    
    % training data
    PosData = Data(:, 1:2);
    VelData = zeros(size(PosData));
    AccData = zeros(size(PosData));
    
    for i = 1:size(PosData, 2)
        VelData(:, i) = gradient(PosData(:, i), dt);
        AccData(:, i) = gradient(VelData(:, i), dt);
    end
    
    PartialData = [PosData, VelData];
    PartialdData = [VelData, AccData];
    
    % test data
    PosData = Data_test(:, 1:2);
    VelData = zeros(size(PosData));
    AccData = zeros(size(PosData));
    
    for i = 1:size(PosData, 2)
        VelData(:, i) = gradient(PosData(:, i), dt);
        AccData(:, i) = gradient(VelData(:, i), dt);
    end
    % VelData = dData_test(:,1:2);
    % AccData = dData_test(:,3:4);
    
    PartialData_test = [PosData, VelData];
    PartialdData_test = [VelData, AccData];
    % 
    % % --- Assign final variables ---
    Data = PartialData;
    dData = PartialdData;
    Data_test = PartialData_test;
    dData_test = PartialdData_test;
    % 
end

if Partial==1 && Smoothing==1
% Parameters for smoothing (tweak as needed)
sgolayOrder  = sgolay_order;   % Polynomial order
sgolayWindow = sgolay_window;

% Preallocate
PosData = Data(:, 1:2);
VelData = zeros(size(PosData));
AccData = zeros(size(PosData));

% Loop over each dimension (e.g., x, y)
for i = 1:size(PosData, 2)
    % Smooth position
    PosData(:, i) = sgolayfilt(PosData(:, i), sgolayOrder, sgolayWindow);

    % Compute smoothed velocity (first derivative)
    VelData(:, i) = gradient(PosData(:, i), dt);

    % Compute smoothed acceleration (second derivative)
    AccData(:, i) = gradient(VelData(:, i), dt);
end

% Collect outputs
PartialData  = [PosData, VelData];
PartialdData = [VelData, AccData];

    % Preallocate
    PosData = Data_test(:, 1:2);

    % Loop over each dimension (e.g., x, y)
    for i = 1:size(PosData, 2)
        % Smooth position
        PosData(:, i) = sgolayfilt(PosData(:, i), sgolayOrder, sgolayWindow);
    
        % Compute smoothed velocity (first derivative)
        VelData(:, i) = gradient(PosData(:, i), dt);
        VelData(:, i) = sgolayfilt(VelData(:, i), sgolayOrder, sgolayWindow);
    
        % Compute smoothed acceleration (second derivative)
        AccData(:, i) = gradient(VelData(:, i), dt);
        AccData(:, i) = sgolayfilt(AccData(:, i), sgolayOrder, sgolayWindow);
    end

    %VelData = dData_test(:,1:2);
    %AccData = dData_test(:,3:4);

    PartialData_test = [PosData, VelData];
    PartialdData_test = [VelData, AccData];

    Data = PartialData;
    dData = PartialdData;
    Data_test = PartialData_test;
    dData_test = PartialdData_test;
end

if Partial==1 && Smoothing==2
    
    % training data
    PosData = Data(:, 1:2);
    % VelData = zeros(size(PosData));
    % AccData = zeros(size(PosData));
    
    for i = 1:size(PosData, 2)
        eyu = TVRegDiff(PosData(:,i), 200, 0.03, [], 'small', 1e-6, dt, 0, 0);
        VelData(:, i) = eyu(1:end-1);
        eyu = TVRegDiff(VelData(:,i), 200, 0.03, [], 'small', 1e-6, dt, 0, 0);
        AccData(:, i) = eyu(1:end-1);
        
    end
    figure(2)
    plot(tspan,AccData(:,1),'linewidth',3,'color','black')
    box('off')
    axis('on')
    set(gca,'FontSize',24)
    
    figure(3)
    plot(tspan,VelData(:,1),'linewidth',3,'color','black')
    box('off')
    axis('on')
    set(gca,'FontSize',24)
    PartialData = [PosData, VelData];
    PartialdData = [VelData, AccData];
    
    % test data
    PosData = Data_test(:, 1:2);
    
    for i = 1:size(PosData, 2)
        VelData(:, i) = gradient(PosData(:, i), dt);
        AccData(:, i) = gradient(VelData(:, i), dt);
    end
     VelData = dData_test(:,1:2);
     AccData = dData_test(:,3:4);
    
    PartialData_test = [PosData, VelData];
    PartialdData_test = [VelData, AccData];
    % 
    % % --- Assign final variables ---
    Data = PartialData;
    dData = PartialdData;
    Data_test = PartialData_test;
    dData_test = PartialdData_test;
    % 
end

    %% === White box regression ===
    % (KEEP your Nonlinlsq block as-is)


    %% === Compare ===
    % (KEEP Compare block with simulation and overlay plots)

    %% === Sindy regression ===
if Sindy == 1

% Get the number of states we have
[dtat_length,n_state]=size(Data);

% Define the control input(Should be zero in our example)
n_control=2;

% Choose whether you want to display actual ODE or not
disp_actual_ode=1;

% If the ODEs you want to display is the actual underlyting dynamics of the
% system, please set actual as 1
actual=1;

% Print the actual ODE we try to discover
% if Control==1
% Print_ODEs(@(t,y,inp)DouPenODE(t, y, inp, l1, l2, m1, m2, Jm1, Jm2, b1, b2, tau1, tau2, tanh_k),n_state,n_control,disp_actual_ode,actual);
% else
% n_control=0;
% Print_ODEs(@(t,y)DouPenODE(t, y, inp, l1, l2, m1, m2, Jm1, Jm2, b1, b2, tau1, tau2, tanh_k),n_state,n_control,disp_actual_ode,actual);
% end

% Create symbolic states
dz=sym('dz',[n_state,1]);

% Now we first create the parameters of the function right hand side
Highest_Poly_Order_Guess=1;
Highest_Trig_Order_Guess=1;
Highest_U_Order_Guess=0;

% Then create the right hand side library parameters
Highest_Poly_Order=1;
Highest_Trig_Order=4;
Highest_U_Order=1;
Highest_dPoly_Order=1;

%% Define parameters for the sparese regression
lam=[1e-4;5e-4;1e-3;2e-3;3e-3;4e-3;5e-3;6e-3;7e-3;8e-3;9e-3;1e-2;2e-2;3e-2;4e-2;5e-2;...
    6e-2;7e-2;8e-2;9e-2;1e-1;2e-1;3e-1;4e-1;5e-1;6e-1;7e-1;8e-1;9e-1;1;1.5;2;2.5;3;3.5;4;4.5;5;...
    6;7;8;9;10;20;30;40;50;100;200];

N_iter=5;
disp=0;
NormalizeLib=1;

for iter=1:n_state
    %fprintf('\n \n Calculating the %i expression...\n',iter)
    
    % According to the previous parameter generate the left hand side guess
    [LHS_Data,LHS_Sym]=GuessLib(Data,dData(:,iter),iter,u,Highest_Poly_Order_Guess,Highest_Trig_Order_Guess,Highest_U_Order_Guess, Integral, dt);
    
    %Generate the corresponding data
    [SINDy_Data,SINDy_Struct]=SINDyLib(Data,dData(:,iter),iter,u,Highest_Poly_Order,Highest_Trig_Order,Highest_U_Order,Highest_dPoly_Order, Integral, dt);
    %input("press space to continue, check LHS and RHS data")
    for i = 1:length(SINDy_Struct)
        %fprintf('%s\n', SINDy_Struct{i});
    end

    % Run the for loop and try all the left hand guess
    for i=1:length(LHS_Sym)
        if iter==1 && i==1
            Xi=cell(n_state,length(LHS_Sym),length(lam));
            ODE=cell(n_state,length(LHS_Sym),length(lam));
            ODEs=cell(n_state,length(LHS_Sym),length(lam));
            %fprintf("xd")
        end

        % Print the left hand side that we are testing
        %fprintf('\t Testing the left hand side as %s:\n',char(LHS_Sym{i}))
        
        % Exclude the guess from SINDy library
        [RHS_Data,RHS_Struct]=ExcludeGuess(SINDy_Data,SINDy_Struct,LHS_Sym{i});
        
        parfor j=1:length(lam)
            % Select the sparse threashold
            lambda=lam(j);
    
            % Perform the sparse regression problem
            [Xi{iter,i,j},ODE{iter,i,j}]=sparsifyDynamics(RHS_Data,LHS_Data(:,i),LHS_Sym{i},lambda,N_iter,RHS_Struct,disp,NormalizeLib);
            
            % Perform sybolic calculation and solve for dX
            digits(6)
            ODE_Guess=vpa(solve(LHS_Sym{i}==ODE{iter,i,j},dz(iter)));
            
            % Print the discovered ODE
            %fprintf(strcat('\t The corresponding ODE we found is: ',char(dz(iter,1)),'=',char((ODE_Guess)),'\n \n'));
            
            % Store the result
            ODEs{iter,i,j}=ODE_Guess;
        end
    end
end

%% Now generate the ODE function file and test the accuracy of the
% identified system

%fprintf('\v Start calculating the best model that could represent the training data...\n \n')
for iter=1:n_state
    % Print which expression are you working on
    %fprintf('\t Calculating the best model for the %d expression...\n',iter)
    
    for i=1:length(LHS_Sym)
        % Print the process
        %fprintf('\t Calculating the score of previously found ODE on the test data, %d %% finished. \n',round((i/length(LHS_Sym))*100))
        
        for j=1:length(lam)
            % If the previous ODE is 0, set the score as NaN, else calculate
            % it.
            if isempty(ODEs{iter,i,j})
                ODE_Not_Exist=1;
                Score(iter,i,j)=NaN;
            else
                % Generate the ODE file
                Generate_ODE_RHS(ODEs{iter,i,j},n_state,n_control);
                % Calculate the accuracy of the file
                Score(iter,i,j)=Get_Score(dData_test(:,iter),Data_test,u_test,Control,tspan,state0_test,Shuffle,iter);
            end
        end
        
        % Get the best lambda
        [minVal1(iter,i),minIndex1(iter,i)]=min(Score(iter,i,:));
        
    end
    
    % Get the best score and use this ODE file
    [minVal2(iter,1),minIndex2(iter,1)]=min(minVal1(iter,:));
    
    % Store the best ODE
    ODE_Best(iter,1)=ODEs{iter,minIndex2(iter,1),minIndex1(iter,minIndex2(iter,1))};
    
    % Print the Result
    %fprintf('\n\n\n\t The SINDy-PI discovered Best ODE for the %d expression is:\n',iter)
    %fprintf('\t %s = %s \n\n\n',char(dz(iter)),char((ODE_Best(iter,1)))')
end


%% Get the simulation result
% Generate the ODE file
%fprintf('\n\n\n\v Generating the Best Model for comparision...\n')
Generate_ODE_RHS(ODE_Best(:,1),n_state,n_control);

% Define a new test data for the comparison
Noise_test = 0;
[dData_Es,Data_Es]=Get_Sim_Data(@(t,z,inp)Sindy_ODE_RHS(t,z,inp),state0_test,u_test,tspan_test,Noise_test,Control,Shuffle);
%% Compute RMSE per state
errors = Data_Es - Data_test;
RMSE = sqrt(mean(errors.^2, 1));  % 1x4 vector, one RMSE per state

% Print nicely
% fprintf('\n=== RMSE Results ===\n');
% fprintf('Theta1:     %.6f\n', RMSE(1));
% fprintf('Theta2:     %.6f\n', RMSE(2));
% fprintf('dTheta1:    %.6f\n', RMSE(3));
% fprintf('dTheta2:    %.6f\n', RMSE(4));

results.RMSE = RMSE;
end