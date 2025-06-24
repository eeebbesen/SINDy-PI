y = [1; 2; 3; 4];
dt = 0.1;

cumtrapz(y) * dt        % WRONG way
cumtrapz(dt, y)         % CORRECT way
dt = 0.001

sgolayWindow = round((0.002 / dt) * 185) + mod(round((0.002 / dt) * 185) + 1, 2)