function [X_est, Z_filt] = runEKF(Z, T, U)
% runEKF  Run EKF filtering on 2xN measurement matrix
%
%   [X_est, Z_filt] = runEKF(Z, T, U)
%
%   Inputs:
%       Z : 2xN matrix of noisy measurements (angles)
%       T : 1xN time vector (optional, default = 1:N)
%       U : control inputs over N steps (optional, default = zeros)
%
%   Outputs:
%       X_est : 4xN filtered state estimates
%       Z_filt: 2xN filtered measurement estimates (angles)

    if nargin < 2 || isempty(T)
        T = 1:size(Z,2); % default time vector
    end
    if nargin < 3
        U = []; % default no inputs
    end

    % --- Load EKF functions ---
    ekfFuncs = EKF_precompute(1); % 1 = regular, ~1 = sindy
    f_k_func = ekfFuncs.f_k;
    h_k_func = ekfFuncs.h_k;
    F_k_func = ekfFuncs.F_k;
    H_k_func = ekfFuncs.H_k;

    % --- Initial Conditions ---
    x_hat_k_plus = sim.initial_state;                 % initial state
    P_k_plus     = diag([0.01, 0.01, 0.01, 0.01]);    % initial covariance
    u_k_prev     = sim.initial_input;                 % prev input
    t            = T(1);

    % Noise covariances (tune!)
    Q = 1e-4 * eye(4);   % process noise
    R = 1e-2 * eye(2);   % measurement noise

    % Allocate storage
    N = size(Z,2);
    X_est = zeros(4,N); % state estimates
    Z_filt = zeros(2,N); % measurement estimates

    % --- Main EKF loop ---
    for k = 1:N
        % Current measurement and input
        z_k = Z(:,k);
        if ~isempty(U)
            u = U(:,k);
        else
            u = zeros(size(u_k_prev));
        end
        time = T(k);

        % Prediction
        delta_t = time - t;
        F_k = double(F_k_func(x_hat_k_plus, u_k_prev, delta_t));
        x_hat_k_minus = double(f_k_func(x_hat_k_plus, u_k_prev, delta_t));
        P_k_minus = F_k * P_k_plus * F_k.' + Q;

        % Update
        H_k = double(H_k_func(x_hat_k_minus));
        K_k = P_k_minus * H_k.' / (H_k * P_k_minus * H_k.' + R);
        x_hat_k_plus = x_hat_k_minus + K_k * (double(z_k(:)) - double(h_k_func(x_hat_k_minus)));
        P_k_plus = (eye(4) - K_k * H_k) * P_k_minus;

        % Save
        X_est(:,k) = x_hat_k_plus;
        Z_filt(:,k) = double(h_k_func(x_hat_k_plus)); % filtered measurements

        % Prepare for next iteration
        t = time;
        u_k_prev = u;
    end
end
