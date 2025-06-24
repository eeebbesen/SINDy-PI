function res = residuals(state_data, u_data, time, l1, l2, tanh_k, m1, m2, b1, b2, tau1, tau2)

    % % Extract measured states
    % theta_measured = state_data(:,1:2);        % Measured θ1, θ2
    % dtheta_measured = state_data(:,3:4);       % Measured dθ1, dθ2
    % 
    % % Initial state vector [θ1, θ2, dθ1, dθ2]
    % state0 = [theta_measured(1,:), dtheta_measured(1,:)];
    % 
    % % Preallocate arrays for simulated states and derivatives
    % y_list = zeros(length(time), length(state0));
    % d_y_list = zeros(length(time), length(state0));
    % 
    % y_list(1,:) = state0;
    % 
    % % Define options for ode15s if needed
    % opts = [];
    % 
    % % Define ODE function handle
    % % Note: u_interp(t) replaced by u(:,i) in the loop to match your method
    % ODE = @(t, y, u) DouPenODE(t, y, u, l1, l2, m1, m2, b1, b2, tau1p, tau1n, tau2p, tau2n, tanh_k);
    % 
    % for i = 2:length(time)
    %     tspan = time(i-1:i);  % Current small time interval
    %     if mod(i, 10000) == 0
    %         disp(i)
    %     end
    %     % Solve ODE over this small interval using previous state as initial condition
    %     try
    %         [~, y_sol] = ode15s(@(t, y) ODE(t, y, u_data(:, i-1)), tspan, state0.', opts);
    %         y_list(i,:) = y_sol(end, :);
    %     catch
    %         y_list(i,:) = Inf(1, length(state0));
    %     end
    % 
    %     % Compute derivative at the new state using current input
    %     d_y_list(i,:) = ODE(0, y_list(i,:).', u_data(:, i));
    % 
    %     % Update initial state for next iteration
    %     state0 = y_list(i,:);
    % end


    noise = 0;
    Control = 1;
    Shuffle = 0;
    state0 = state_data(1,1:4).';
    [dData,Data] = Get_Sim_Data(@(t,y,inp)DouPenODE(t, y, inp, l1, l2, m1, m2, b1, b2, tau1, tau2, tanh_k),state0,u_data,time,noise,Control,Shuffle);

    % Residuals between simulated and measured θ1, θ2 (flattened for lsqnonlin)
    res(:,1) = (Data(:,1) - state_data(:,1));
    res(:,2) = (Data(:,2) - state_data(:,2));
    

    % estimated_data1 = DouPenODE(0, state_data(:,1:4).', u_data, l1, l2, m1, m2, b1, b2, tau1p, tau1n, tau2p, tau2n, tanh_k);
    % estimated_data = estimated_data1.';
    % residuals = (estimated_data(:,1) - state_data(:,3)).^2 + (estimated_data(:,2) - state_data(:,4)).^2 + ...
    %             (estimated_data(:,3) - state_data(:,5)).^2 + (estimated_data(:,4) - state_data(:,6)).^2;
    % res = sum(residuals);

end
