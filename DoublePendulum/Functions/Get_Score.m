%% This function will calculate the l2 norm error between the discovered system and actual test data
% Last Updated: 2019/04/22
% Coded By: K

function [Score]=Get_Score(dData_test,Data_test,u,Control,tspan,state0,Shuffle,iter)
% % Get the ODE simulation result
% Noise=0;
% 
% % Using ODEs to get the simulation data
% [dData_Es,Data_Es]=Get_Sim_Data(@(t,z,u)Sindy_ODE_RHS(t,z,u),state0,u,tspan,Noise,Control,Shuffle);
% 
% % Get the score of the result
% [n,m]=size(Data_Es);
% Score=sum(norm(dData_test(1:n,:)-dData_Es,"fro"))+sum(norm(Data_test(1:n,:)-Data_Es,"fro"));


%% Using one step prediction to get the simulation data
ind = iter;
dt = tspan(1,2)-tspan(1,1);
dData_Es = zeros(size(dData_test));
for i=1:length(Data_test(:,ind))
dData_Es(i)=Sindy_ODE_RHS(0,Data_test(i,:)',u(:,i));
Data_test(i+1,ind) = Data_test(i,ind) + dData_Es(i)*dt;
end

%dData_Es=Sindy_ODE_RHS(0,Data_test',u)';

% Get the score of the result
Score=sum(norm(dData_test-dData_Es));





