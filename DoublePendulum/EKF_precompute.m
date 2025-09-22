function ekfFuncs = EKF_precompute(choice)
% Precomputes symbolic functions and Jacobians for EKF from saved symbolic ODE
% Inputs:
%   - odeMatFile: string path to .mat file with symbolic ODE (expects variable 'f_cont_sym')
%   - params: struct with at least l1 and l2
% Output:
%   - ekfFuncs: struct with EKF function handles + symbolic forms
choice = 1;
    % === 1. Load symbolic ODE ===
    if choice == 1
    data = load("saved_regode.mat");
    else
    data = load("saved_sindyode.mat"); 
    end
    f_cont = data.f;

    % === 2. Declare symbolic variables ===
    syms theta1 theta2 dtheta1 dtheta2 u1 u2 dt real
    x = [theta1; theta2; dtheta1; dtheta2];
    u = [u1; u2];
    
    if choice ~=1
    % Define substitution mapping
    syms z1 z2 z3 z4 real
    z = [z1; z2; z3; z4];
    
    % Replace z1–z4 with actual state variables
    f_cont = subs(f_cont, z, x);
    end

    % === 3. Euler Discretization ===
    f_k_sym = x + dt * f_cont;

    % === 4. Measurement model ===

    h_k_sym = [
        theta1;
        theta2;
    ];

    % === 5. Jacobians ===
    F_k_sym = jacobian(f_k_sym, x);
    H_k_sym = jacobian(h_k_sym, x);

    % === 6. Convert to MATLAB functions ===
    ekfFuncs.f_k = matlabFunction(f_k_sym, 'Vars', {x, u, dt}, 'File', 'f_k_func');
    ekfFuncs.h_k = matlabFunction(h_k_sym, 'Vars', {x}, 'File', 'h_k_func');
    ekfFuncs.F_k = matlabFunction(F_k_sym, 'Vars', {x, u, dt}, 'File', 'F_k_func');
    ekfFuncs.H_k = matlabFunction(H_k_sym, 'Vars', {x}, 'File', 'H_k_func');

    % Optional: include symbolic forms
    ekfFuncs.f_k_sym = f_k_sym;
    ekfFuncs.h_k_sym = h_k_sym;
    ekfFuncs.F_k_sym = F_k_sym;
    ekfFuncs.H_k_sym = H_k_sym;
end
