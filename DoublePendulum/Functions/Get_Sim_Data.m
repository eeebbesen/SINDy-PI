function [d_Data,Data]=Get_Sim_Data(ODE,state0,u,tspan,Noise,Control,Shuffle)
%% Get the size of the state and control
[N1,M1]=size(state0);
[N2,M2]=size(u);

%% ODE options
opts = odeset('RelTol', 1e-3, 'AbsTol', 1e-6, 'MaxStep', 0.01);  % limit max step for stiff regions

%% Determine the left hand side derivative
if Control==1
    y_list(1,:)=state0;
    d_y_list(1,:)=ODE(0,y_list(1,:).',u(:,1));

    for i=2:length(u)
         if mod(i, 10) == 0
             %disp(['Step ' num2str(i) ' / ' num2str(length(u))]);
         end
         
         % Step-by-step integration for discrete input with freeze handling
         u_step = u(:, i-1);
         ode_for_solver = @(t,y) ODE_wrapper(t, y, ODE, u_step);

         [~, y_1] = ode113(ode_for_solver, tspan(i-1:i), y_list(i-1,:).', opts);

         y_list(i,:) = y_1(end,:);
         d_y_list(i,:) = ODE_wrapper(0, y_list(i,:).', ODE, u(:,i));
    end
else
    % Simulate the system ODE
    [t,y]=ode45(@(t,y)ODE(t,y),tspan,state0);
    y_list=y;

    dt = tspan(2) - tspan(1);
    VelData = zeros(size(y_list));
    for i = 1:2
        VelData(:, i) = gradient(y(:, i), dt);
    end

    % Get the derivative data
    d_y_list=ODE(0,y_list')';
    y_list = [y(:,1),y(:,2),VelData(:,3),VelData(:,4)];
end

%% Add some noise to the system
for i=1:N1
    Data(:,i)=y_list(:,i)+Noise*randn(size(y_list(:,i)));
end

for i=1:N1
    d_Data(:,i)=d_y_list(:,i)+Noise*randn(size(d_y_list(:,i)));
end

%% Shuffle the data
if Shuffle==1
    Sequence=randperm(size(Data,1));
    Data=Data(Sequence,:);
    d_Data=d_Data(Sequence,:);
end

end

%% --- Helper wrapper for freezing low-velocity / low-control states ---
function dy = ODE_wrapper(~, y, ODE, u_step)
    dy = ODE(0, y, u_step);
    
    vel_idx = 3:4; % velocity components
    if norm(y(vel_idx)) < 1e-8 && norm(u_step) < 1e-8
        dy = zeros(size(dy)); % freeze near-rest
    end
    
    % Optional clamp to prevent extreme spikes
    dy = max(min(dy, 1e4), -1e4);
end
