%% EKF
ekfFuncs = EKF_precompute(1); % 1 for reg, ~1 for sindy
 % self.x_hat_k_plus = sim.initial_state
 %        self.P_k_plus = np.diag([0.01,0.01,0.01,0.01])
 %        self.x_hat_k_prev_plus = 0
 %        self.P_k_prev_plus = 0
 %        self.u_k_prev = sim.initial_input

f_k_func = ekfFuncs.f_k;
h_k_func = ekfFuncs.h_k;
F_K_func = ekfFuncs.F_k;
H_k_func = ekfFuncs.H_k;

% ChattGPTAD kod vvv
% Update previous estimates
x_hat_k_prev_plus = x_hat_k_plus;
P_k_prev_plus = P_k_plus;
k_count = k_count + 1;

delta_t = time - t;

% Calculate F_k (Jacobian of the state transition function)
F_k = double(F_k_func(x_hat_k_prev_plus, u_k_prev, delta_t));  % Assuming F_k_func is a handle

% Predict state using nonlinear function f_k
x_hat_k_minus = double(f_k_func(x_hat_k_prev_plus, u_k_prev, delta_t));  % Flatten not needed

% Calculate P_k_minus
P_k_minus = F_k * P_k_prev_plus * F_k.' + Q;

% Calculate H_k (Jacobian of the measurement function)
H_k = double(H_k_func(x_hat_k_minus));  % Assuming H_k_func is a handle

% Calculate Kalman Gain
K_k = P_k_minus * H_k.' / (H_k * P_k_minus * H_k.' + R);

% Update state estimate using measurement z_k
x_hat_k_plus = x_hat_k_minus + K_k * (double(z_k(:)) - double(h_k_func(x_hat_k_minus)));

% Update error covariance
P_k_plus = (eye(4) - K_k * H_k) * P_k_minus;  % Identity matrix of size 4

% Save current values for next iteration
t = time;
u_k_prev = u;
