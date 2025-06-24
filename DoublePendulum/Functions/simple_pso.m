function [best_pos, best_val] = simple_pso(objfun, lb, ub, nVars, options)
% SIMPLE_PSO A basic Particle Swarm Optimizer for bound-constrained problems
%
% Inputs:
%   objfun - function handle returning scalar cost
%   lb, ub - lower and upper bounds (vectors)
%   nVars  - number of variables
%   options - struct with fields:
%       .MaxIter (max iterations)
%       .SwarmSize (number of particles)
%       .Inertia, .Cognitive, .Social (PSO parameters)
%
% Outputs:
%   best_pos - best found position
%   best_val - objective value at best_pos

% Default options if not provided
if nargin < 5
    options.MaxIter = 200;
    options.SwarmSize = 50;
    options.Inertia = 0.7;
    options.Cognitive = 1.4;
    options.Social = 1.4;
end

% Initialize particles
pos = repmat(lb, options.SwarmSize, 1) + rand(options.SwarmSize, nVars) .* (repmat(ub - lb, options.SwarmSize,1));
vel = zeros(options.SwarmSize, nVars);

% Evaluate initial particles
fitness = zeros(options.SwarmSize, 1);
for i = 1:options.SwarmSize
    fitness(i) = objfun(pos(i,:));
end

% Personal bests
pbest_pos = pos;
pbest_val = fitness;

% Global best
[best_val, idx] = min(pbest_val);
best_pos = pbest_pos(idx,:);
disp("hello")
% Main loop
for iter = 1:options.MaxIter
    for i = 1:options.SwarmSize
        % Update velocity
        r1 = rand(1, nVars);
        r2 = rand(1, nVars);
        vel(i,:) = options.Inertia * vel(i,:) ...
            + options.Cognitive * r1 .* (pbest_pos(i,:) - pos(i,:)) ...
            + options.Social * r2 .* (best_pos - pos(i,:));
        
        % Update position
        pos(i,:) = pos(i,:) + vel(i,:);
        
        % Enforce bounds
        pos(i,:) = max(pos(i,:), lb);
        pos(i,:) = min(pos(i,:), ub);
        
        % Evaluate
        fitness(i) = objfun(pos(i,:));
        
        % Update personal best
        if fitness(i) < pbest_val(i)
            pbest_val(i) = fitness(i);
            pbest_pos(i,:) = pos(i,:);
        end
    end
    
    % Update global best
    [current_best_val, idx] = min(pbest_val);
    if current_best_val < best_val
        best_val = current_best_val;
        best_pos = pbest_pos(idx,:);
    end
    
    % Display progress
    fprintf('Iter %3d | Best Cost: %.6f\n', iter, best_val);
end

end
