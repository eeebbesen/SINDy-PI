function dydt = DouPenODE(t, y, u, l1, l2, m1, m2, Jm1, Jm2, b1, b2, tau1, tau2, tanh_k, usesym)
    % theta1 = y(1,:);
    % theta2 = y(2,:);
    % dtheta1 = y(3,:);
    % dtheta2 = y(4,:);
    % u1 = u(1,:);
    % u2 = u(2,:);
    % 
    % % Precompute terms for readability
    % cos_theta2 = cos(theta2);
    % sin_theta2 = sin(theta2);
    % tanh_dtheta1 = tau1*tanh(tanh_k * dtheta1);
    % tanh_dtheta2 = tau2*tanh(tanh_k * dtheta2);
    % 
    % % Denominators
    %         % (Jm1 * Jm2 + Jm1 * l2**2 * m2 + Jm2 * l1**2 * m1 +
    %         %  Jm2 * l1**2 * m2 + 2 * Jm2 * l1 * l2 * m2 * cos(theta2) +
    %         %  Jm2 * l2**2 * m2 + l1**2 * l2**2 * m1 * m2 -
    %         %  l1**2 * l2**2 * m2**2 * cos(theta2)**2 + l1**2 * l2**2 * m2**2)
    % 
    % D1 = l1^2 * l2 * m1 - l1^2 * l2 * m2 * cos_theta2.^2 + l1^2 * l2 * m2;
    % %D1 = Jm1*Jm2 + Jm1*l2^2*m2 + Jm2*l1^2*m1 + Jm2*l1^2*m2 + 2*Jm2*l2*l2*m2*cos_theta2 + Jm2*l2^2*m2 + l1^2*l2^2*m1*m2 - l1^2*l2^2*m2^2*cos_theta2.^2 + l1^2*l2^2*m2^2;
    % % disp(D1)
    % % disp(l1^2)
    % % disp(l2)
    % % disp(m1)
    % % disp(m2)
    % % disp(cos_theta2.^2)
    % 
    % %         (Jm1 * Jm2 + Jm1 * l2**2 * m2 + Jm2 * l1**2 * m1 +
    % %          Jm2 * l1**2 * m2 + 2 * Jm2 * l1 * l2 * m2 * cos(theta2) +
    % %          Jm2 * l2**2 * m2 + l1**2 * l2**2 * m1 * m2 -
    % %          l1**2 * l2**2 * m2**2 * cos(theta2)**2 + l1**2 * l2**2 * m2**2)
    % % 
    % D2 = l1^2 * l2^2 * m1 * m2 - l1^2 * l2^2 * m2.^2 .* cos_theta2.^2 + l1^2 * l2^2 * m2.^2;
    % %D2 = Jm1*Jm2 + Jm1*l2^2*m2 + Jm2*l1^2*m1 + Jm2*l1^2*m2 + 2*Jm2*l1*l2*m2*cos_theta2 + Jm2*l2^2*m2 + l1^2*l2^2*m1*m2 - l1^2*l2^2*m2^2*cos_theta2.^2 + l1^2*l2^2*m2^2;
    % 
    % % ddtheta1
    % num1 = (-l1 * cos_theta2 - l2) .* (-b2 .* dtheta2 - dtheta1.^2 .* l1 .* l2 .* m2 .* sin_theta2 - tanh_dtheta2 + u2);
    % num2 = l2 .* (-b1 .* dtheta1 + 2 .* dtheta1 .* dtheta2 .* l1 .* l2 .* m2 .* sin_theta2 + ...
    %        dtheta2.^2 .* l1 .* l2 .* m2 .* sin_theta2 - tanh_dtheta1 + u1);
    % % 
    % % num1 = (Jm2 + l2**2 * m2) * (-b1 * dtheta1 +
    % %                              2 * dtheta1 * dtheta2 * l1 * l2 * m2 * sin(theta2) +
    % %                              dtheta2**2 * l1 * l2 * m2 * sin(theta2) -
    % %                              tau1 * tanh(tanh_k * dtheta1) + u1)
    % % 
    % % num2 = (-l1 * l2 * m2 * cos(theta2) - l2**2 * m2) * (
    % %     -b2 * dtheta2 - dtheta1**2 * l1 * l2 * m2 * sin(theta2) -
    % %     tau2 * tanh(tanh_k * dtheta2) + u2)
    % 
    % % num1 = (Jm2 + l2^2*m2).*(-b1*dtheta1 + 2*dtheta1.*dtheta2*l1*l2*m2.*sin_theta2 + dtheta2.^2*l1*l2*m2.*sin_theta2 - tau1*tanh_dtheta1 + u1);
    % % num2 = (-l1*l2*m2*cos_theta2 - l2^2*m2).*(-b2*dtheta2 - dtheta1.^2*l1*l2*m2.*sin_theta2 - tau2*tanh_dtheta2 + u2);
    % ddtheta1 = (num1 + num2)./ D1;
    % 
    % % disp(num1)
    % % disp(num2)
    % 
    % % ddtheta2
    % num2_1 = m2 .* l2 .* ((-l1 .* cos_theta2 - l2) .* ...
    %           (-b1 .* dtheta1 + 2 .* dtheta1 .* dtheta2 .* l1 .* l2 .* m2 .* sin_theta2 + ...
    %            dtheta2.^2 .* l1 .* l2 .* m2 .* sin_theta2 - tanh_dtheta1 + u1));
    % num2_2 = (-b2 .* dtheta2 - dtheta1.^2 .* l1 .* l2 .* m2 .* sin_theta2 - tanh_dtheta2 + u2) .* ...
    %          (l1.^2 .* m1 + l1.^2 .* m2 + 2 .* l1 .* l2 .* m2 .* cos_theta2 + l2.^2 .* m2);
    % 
    % % 
    % % num1 = (-l1 * l2 * m2 * cos(theta2) - l2**2 * m2) * (
    % %     -b1 * dtheta1 + 2 * dtheta1 * dtheta2 * l1 * l2 * m2 * sin(theta2) +
    % %     dtheta2**2 * l1 * l2 * m2 * sin(theta2) -
    % %     tau1 * tanh(tanh_k * dtheta1) + u1)
    % % 
    % % num2 = (-b2 * dtheta2 - dtheta1**2 * l1 * l2 * m2 * sin(theta2) -
    % %         tau2 * tanh(tanh_k * dtheta2) + u2) * (
    % %     Jm1 + l1**2 * m1 + l1**2 * m2 +
    % %     2 * l1 * l2 * m2 * cos(theta2) + l2**2 * m2)
    % 
    % % num2_1 = (-l1*l2*m2*cos_theta2 - l2^2*m2).*(-b1*dtheta1 + 2*dtheta1.*dtheta2*l1*l2*m2.*sin_theta2 + dtheta2.^2*l1*l2*m2.*sin_theta2 - tau1*tanh_dtheta1 + u1);
    % % num2_2 = (-b2*dtheta2 - dtheta1.^2*l1*l2*m2.*sin_theta2 - tau2*tanh_dtheta2 + u2).*(Jm1 + l1^2*m1 + l1^2*m2 + 2*l1*l2*m2*cos_theta2 + l2^2*m2);
    % ddtheta2 = (num2_1 + num2_2) ./ D2;
    % % disp(size(dtheta1))
    % % disp(size(dtheta2))
    % % disp(ddtheta1)
    % % disp(ddtheta2)
    % % Return dydt

% Unpack components
    if usesym == 1
        syms theta1 theta2 dtheta1 dtheta2 u1 u2 real
    else
    theta1 = y(1,:);
    theta2 = y(2,:);
    dtheta1 = y(3,:);
    dtheta2 = y(4,:);
    u1 = u(1,:);
    u2 = u(2,:);
    end
    tanh_dtheta1 = tau1*tanh(tanh_k * dtheta1);
    tanh_dtheta2 = tau2*tanh(tanh_k * dtheta2);
    

% Mass matrix M (element-wise operations)
M11 = (l1.^2 .* (m1 + m2) + l2.^2 .* m2 + 2 .* l1 .* l2 .* m2 .* cos(theta2)) +Jm1;
M12 = m2 .* (l2.^2 + l1 .* l2 .* cos(theta2));
M22 = ones(size(theta1)).*l2.^2 .* m2 + Jm2;

if usesym == 0
N = length(theta1);  % number of samples
end

if usesym == 1
    M = [M11,M12;M12,M22]; % symbolic zeros array
else
    M = zeros(2, 2, N);      % numeric zeros array
    M(1,1,:) = M11;
    M(1,2,:) = M12;
    M(2,1,:) = M12;
    M(2,2,:) = M22;
end



% Coriolis vector C
C1 = -2 .* l1 .* l2 .* m2 .* dtheta1 .* dtheta2 .* sin(theta2) ...
     - l1 .* l2 .* m2 .* dtheta2.^2 .* sin(theta2);
C2 = l1 .* l2 .* m2 .* dtheta1.^2 .* sin(theta2);


if usesym == 1
    C = [C1;C2]; % symbolic zeros array
else
    C = zeros(2, 1, N);      % numeric zeros array
    C(1,1,:) = C1;
    C(2,1,:) = C2;  
end

% Friction vector F
F1 = b1 .* dtheta1 + tanh_dtheta1;
F2 = b2 .* dtheta2 + tanh_dtheta2;

if usesym == 1
    F = [F1;F2]; % symbolic zeros array
else
    F = zeros(2, 1, N);      % numeric zeros array
    F(1,1,:) = F1;
    F(2,1,:) = F2;
end

% Torque vector T
T = [u1; u2];

% Solve for ddtheta (angular accelerations)
if usesym == 1
    ddtheta = M \ (T-C-F);
else
    N = size(y, 2);
    ddtheta = zeros(2, N);
    for i = 1:N
        ddtheta(:,i) = M(:,:,i) \ (T(:,i) - C(:,i) - F(:,i));
    end
end

    dydt = [
        dtheta1;
        dtheta2;
        ddtheta(1,:);
        ddtheta(2,:);
    ];
    % dydt = [
    %     dtheta1;
    %     dtheta2;
    %     ddtheta1;
    %     ddtheta2
    % ];
end


% function dydt=DouPenODE(t,y,u,m1,m2,a1,a2,L1,I1,I2,k1,k2,g)
% dydt=[
%     y(3,:);
%     y(4,:);
%     (-(u(1,:) + I2*k1*y(3,:) + I2*k2*y(3,:) - I2*k2*y(4,:) + a2^2*k1*m2*y(3,:) + a2^2*k2*m2*y(3,:) - a2^2*k2*m2*y(4,:) + L1*a2^3*m2^2*y(4,:).^2.*sin(y(1,:) - y(2,:)) - (L1*a2^2*g*m2^2*sin(y(1,:)))/2 - I2*L1*g*m2*sin(y(1,:)) - (L1*a2^2*g*m2^2*sin(y(1,:) - 2*y(2,:)))/2 - I2*a1*g*m1*sin(y(1,:)) + (L1^2*a2^2*m2^2*y(3,:).^2.*sin(2*y(1,:) - 2*y(2,:)))/2 + L1*a2*k2*m2*y(3,:).*cos(y(1,:) - y(2,:)) - L1*a2*k2*m2*y(4,:).*cos(y(1,:) - y(2,:)) - a1*a2^2*g*m1*m2*sin(y(1,:)) + I2*L1*a2*m2*y(4,:).^2.*sin(y(1,:) - y(2,:)))./(I1*I2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + I2*a1^2*m1 + I1*a2^2*m2 - L1^2*a2^2*m2^2*cos(y(1,:) - y(2,:)).^2 + a1^2*a2^2*m1*m2));
%     (( u(2,:) + I1*k2*y(3,:) + I1*k2*y(4,:) + L1^2*k2*m2*y(3,:) - L1^2*k2*m2*y(4,:) + a1^2*k2*m1*y(3,:) - a1^2*k2*m1*y(4,:) + L1^3*a2*m2^2*y(3,:).^2.*sin(y(1,:) - y(2,:)) + L1^2*a2*g*m2^2*sin(y(2,:)) + I1*a2*g*m2*sin(y(2,:)) + (L1^2*a2^2*m2^2*y(4,:).^2.*sin(2*y(1,:) - 2*y(2,:)))/2 + L1*a2*k1*m2*y(3,:).*cos(y(1,:) - y(2,:)) + L1*a2*k2*m2*y(3,:).*cos(y(1,:) - y(2,:)) - L1*a2*k2*m2*y(4,:).*cos(y(1,:) - y(2,:)) - L1^2*a2*g*m2^2*cos(y(1,:) - y(2,:)).*sin(y(1,:)) + a1^2*a2*g*m1*m2*sin(y(2,:)) + I1*L1*a2*m2*y(3,:).^2.*sin(y(1,:) - y(2,:)) + L1*a1^2*a2*m1*m2*y(3,:).^2.*sin(y(1,:) - y(2,:)) - L1*a1*a2*g*m1*m2*cos(y(1,:) - y(2,:)).*sin(y(1,:)))./(I1*I2 + L1^2*a2^2*m2^2 + I2*L1^2*m2 + I2*a1^2*m1 + I1*a2^2*m2 - L1^2*a2^2*m2^2*cos(y(1,:) - y(2,:)).^2 + a1^2*a2^2*m1*m2))
%     ];
