y = [1; 2; 3; 4];
dt = 0.1;

cumtrapz(y) * dt        % WRONG way
cumtrapz(dt, y)         % CORRECT way