% Modelação e identificação do drone CRAZYFLIE
%% Constantes das Forças e Momentos
%close all
clear all
% initPlots

g = 9.807;
m = 0.035;
Cd = 0.2056;
zI = [0; 0; 1];

% Time vector
dt = 0.01;
t_end = 20;
t = 1:dt:t_end;
NSim = length(t);
tspan = [0 t_end];

% Initial conditions
x0 = [0; 0; 0; 0; 0; 0];
% x_ref = [1;1;1;0;0;0];
x_ref = [1;1;1;0;0;0] * ones(size(t)) - 0.5*[1;1;1;0;0;0]*(t>=10);

%% Linear model

nx = 6;
ny = 6;

A = [zeros(3) eye(3)
     zeros(3) -Cd*eye(3)];
B = [zeros(3)
     eye(3)  ];
C = eye(6);
D = zeros(6,3);

sys = ss(A,B,C,D);


% escalar_p = 0.01;
% escalar_v = 0.01;
% Q = kron(diag([escalar_p, escalar_v]), eye(nx/2));
% escalar_R = 0.01;
% R = escalar_R*diag([1.3759e-06,1.2508e-06,1.6314e-06]); %covariancias obtidas no 1.1 (drone parado) do Lab2

escalar_p = 0.01;
escalar_v = 0.01;
Q = kron(diag([escalar_p, escalar_v]), eye(nx/2));
escalar_R = 0.01;
R = escalar_R*eye(3);

[K,S,P] = lqr(sys,Q,R); % Since N is not specified, lqr sets N to 0.

e = eig(A-B*K);

% Linear simulation
x_lin = zeros(nx,NSim);
x_lin(:,1) = x0;
u = zeros(3,NSim);

for k=1:NSim-1
    y(:,k) = C*x_lin(:,k);
    u(:,k) = -K * (x_lin(:,k)-x_ref(:,k));

    dxdt = A*x_lin(:,k) + B*u(:,k);
    
    % integrate state
    x_lin(:,k+1) = x_lin(:,k) + dt * dxdt;
end


%% Nonlinear Model

x_nlin = zeros(nx, 1);
x_nlin(:,1) = x0;

for k = 1:length(t)-1
    % prepare variables:
    ua = -K*(x_nlin(:,k)-x_ref(:,k));
    
    
    % compute state derivative:
    dpdt = x_nlin(4:6,k);
    dvdt = ua - Cd * x_nlin(4:6,k);

    dxdt = [dpdt; dvdt];
    
    % integrate state
    x_nlin(:,k+1) = x_nlin(:,k) + dt * dxdt;
end

%% Plot (Controlador LQR em Modelo Linear vs Não Linear)

figure;

% Linha 1: p_x e v_x
subplot(3,2,1); hold on; grid on;
plot(t, x_nlin(1,:), 'b'); plot(t, x_lin(1,:), 'r--');plot(t, x_ref(1,:));
ylabel('p_x [m]'); title('Posições');

subplot(3,2,2); hold on; grid on;
plot(t, x_nlin(4,:), 'b'); plot(t, x_lin(4,:), 'r--');plot(t, x_ref(4,:));
ylabel('v_x [m/s]'); title('Velocidades');

% Linha 2: p_y e v_y
subplot(3,2,3); hold on; grid on;
plot(t, x_nlin(2,:), 'b'); plot(t, x_lin(2,:), 'r--');plot(t, x_ref(2,:));
ylabel('p_y [m]');

subplot(3,2,4); hold on; grid on;
plot(t, x_nlin(5,:), 'b'); plot(t, x_lin(5,:), 'r--');plot(t, x_ref(5,:));
ylabel('v_y [m/s]');

% Linha 3: p_z e v_z
subplot(3,2,5); hold on; grid on;
plot(t, x_nlin(3,:), 'b'); plot(t, x_lin(3,:), 'r--');plot(t, x_ref(3,:));
ylabel('p_z [m]'); xlabel('Tempo [s]');

subplot(3,2,6); hold on; grid on;
plot(t, x_nlin(6,:), 'b'); plot(t, x_lin(6,:), 'r--');plot(t, x_ref(6,:));
ylabel('v_z [m/s]'); xlabel('Tempo [s]');

legend('Não Linear','Linear','Referência');
sgtitle('Controlador LQR em Modelo Linear vs Não Linear');


%%  Error Model

x_er = zeros(nx, 1);
x_er(:,1) = x0;


% T = [0;0;0.4];
T=m*g;
Ro = eye(3); 

for k = 1:length(t)-1
    % prepare variables:
    p = x_lin(1:3,k);
    v = x_lin(4:6,k);
    p_d = x_ref(1:3,k);
    v_d = x_ref(4:6,k);

    % ep = p-pd  % ev = v - v_d
    x_er(:,k) = x_lin(:,k) - x_ref(:,k);
 
    % compute state derivative:
    depdt = x_er(4:6,k);
    devdt = -(g*zI + Cd*v) + T/m * zI - v_d;

    dxdt = [depdt; devdt];
    
    % integrate state
    x_er(:,k+1) = x_er(:,k) + dt * dxdt;
end

%% Plot (Controlador LQR em Modelo Linear vs de Erro)

figure;

% Linha 1: p_x e v_x
subplot(3,2,1); hold on; grid on;
plot(t, x_er(1,:), 'b'); plot(t, x_lin(1,:), 'r--');plot(t, x_ref(1,:));
ylabel('p_x [m]'); title('Posições');

subplot(3,2,2); hold on; grid on;
plot(t, x_er(4,:), 'b'); plot(t, x_lin(4,:), 'r--');plot(t, x_ref(4,:));
ylabel('v_x [m/s]'); title('Velocidades');

% Linha 2: p_y e v_y
subplot(3,2,3); hold on; grid on;
plot(t, x_er(2,:), 'b'); plot(t, x_lin(2,:), 'r--');plot(t, x_ref(2,:));
ylabel('p_y [m]');

subplot(3,2,4); hold on; grid on;
plot(t, x_er(5,:), 'b'); plot(t, x_lin(5,:), 'r--');plot(t, x_ref(5,:));
ylabel('v_y [m/s]');

% Linha 3: p_z e v_z
subplot(3,2,5); hold on; grid on;
plot(t, x_er(3,:), 'b'); plot(t, x_lin(3,:), 'r--');plot(t, x_ref(3,:));
ylabel('p_z [m]'); xlabel('Tempo [s]');

subplot(3,2,6); hold on; grid on;
plot(t, x_er(6,:), 'b'); plot(t, x_lin(6,:), 'r--');plot(t, x_ref(6,:));
ylabel('v_z [m/s]'); xlabel('Tempo [s]');

legend('Modelo de Erro','Modelo Linear','Referência');
sgtitle('Controlador LQR em Modelo Linear vs de Erro');

%% Nonlinear Control and Trials
% Lyapunov theory

kp_values = [0.5, 1.0, 2.0 2.5];
kv_values = [0.5 1.0, 1.5, 2.5];
results = struct;
colors = lines(length(kp_values)*length(kv_values));

figure;
subplot(1,1,1); hold on; grid on;

idx = 1;

for i = 1:length(kp_values)
    for j = 1:length(kv_values)

        kp = kp_values(i);
        kv = kv_values(j);

        x_lyap = zeros(nx, NSim);
        x_lyap(:,1) = x0;

        for k = 1:NSim-1
            % prepare variables:
            p = x_lyap(1:3,k);
            v = x_lyap(4:6,k);
            p_d = x_ref(1:3,k);
            v_d = x_ref(4:6,k);
        
            ep = p - p_d;
            ev = v - v_d;
        
            % control law
            u = - kp*ep - kv*ev;
            ua = u + Cd*v + g*zI + v_d;
        
            % compute state derivative:
            dpdt = v;
            dvdt = ua - Cd*v - g*zI;
        
            dxdt = [dpdt; dvdt];
        
            % integrate state
            x_lyap(:,k+1) = x_lyap(:,k) + dt * dxdt;
        end

        % Plot p_x
        plot(t, x_lyap(1,:), 'Color', colors(idx,:), 'DisplayName', ['K_p=', num2str(kp), ', K_v=', num2str(kv)]);
        idx = idx + 1;

        key = sprintf('kp_%d_kv_%d', kp*10, kv*10);
        results.(key).kp = kp;
        results.(key).kv = kv;
        results.(key).x_lyap = x_lyap;
    end
end

% Referência
plot(t, x_ref(1,:), 'k--', 'DisplayName', 'Referência');

xlabel('Tempo [s]');
ylabel('p_x [m]');
title('Evolução de p_x para diferentes (K_p, K_v)');
legend('Location', 'best');


%% plot (Controlador Não Linear (Lyapunov))

kp_sel = 2.5;
kv_sel = 2.5;

key_sel = sprintf('kp_%d_kv_%d', kp_sel*10, kv_sel*10);
x_lyap = results.(key_sel).x_lyap;



figure;

% Linha 1: p_x e v_x
subplot(3,2,1); hold on; grid on;
plot(t, x_lyap(1,:), 'b');plot(t, x_ref(1,:));
ylabel('p_x [m]'); title('Posições');

subplot(3,2,2); hold on; grid on;
plot(t, x_lyap(4,:), 'b');plot(t, x_ref(4,:));
ylabel('v_x [m/s]'); title('Velocidades');

% Linha 2: p_y e v_y
subplot(3,2,3); hold on; grid on;
plot(t, x_lyap(2,:), 'b');plot(t, x_ref(2,:));
ylabel('p_y [m]');

subplot(3,2,4); hold on; grid on;
plot(t, x_lyap(5,:), 'b');plot(t, x_ref(5,:));
ylabel('v_y [m/s]');

% Linha 3: p_z e v_z
subplot(3,2,5); hold on; grid on;
plot(t, x_lyap(3,:), 'b');plot(t, x_ref(3,:));
ylabel('p_z [m]'); xlabel('Tempo [s]');

subplot(3,2,6); hold on; grid on;
plot(t, x_lyap(6,:), 'b');plot(t, x_ref(6,:));
ylabel('v_z [m/s]'); xlabel('Tempo [s]');
legend('(Lyapunov)', 'Referência');
sgtitle('Controlador Não Linear obtido pela análise de Lyapunov');





%% Motion Planning
% RTT*

% Define flying arena 
PFly = [-1.2 1.2; -2.1 2.1; 0.0 2.0];

% Define obstacles
PNoFly1 = [-0.6,0.7,0.0,1.2,0.1,1.2];
PNoFly2 = [-0.6,-0.8,0.6,1.2,0.1,1.4];

% Define initial positions
P_Start = [-0.3 0.3; -0.3 0.3; 0 0];
PStart = randInBounds(P_Start);

% Define goal positions
P_Goal1 = [-0.2 0.2; 0.9 1.3; 0.4 0.8];
PGoal1 = randInBounds(P_Goal1);
P_Goal2 = [-0.2 0.2; -1.3 -0.9; 1.0 1.4];
PGoal2 = randInBounds(P_Goal2);

EPS = 0.2;
numNodes = 2000;


% plot
figure
view(3);
axis([PFly(1,:) PFly(2,:) PFly(3,:)])
hold on
patch(obsgen(PNoFly1))
patch(obsgen(PNoFly2))

% Plot dos pontos
plot3(PStart(1), PStart(2), PStart(3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g')
plot3(PGoal1(1), PGoal1(2), PGoal1(3), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm')
plot3(PGoal2(1), PGoal2(2), PGoal2(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')

[traj1, nodes1] = RRTStar_3D(PStart, PGoal1, PNoFly1, PNoFly2, EPS, numNodes, PFly);
[traj2, nodes2] = RRTStar_3D(PGoal1, PGoal2, PNoFly1, PNoFly2, EPS, numNodes, PFly);

    
plotPath(traj1,'b');  % Start → Goal1
plotPath(traj2,'b');  % Goal1 → Goal2
title('Trajetória completa: Start → Goal1 → Goal2')





%%


trajTotal = [traj1; traj2]';
NSim2 = size(trajTotal, 2);
trajTotal = [trajTotal;zeros(3,NSim2)];

t_end2 = 50;

dt2 = t_end2 / NSim2;
t2 = 1:dt2:t_end2;

% Linear simulation
x = zeros(nx,NSim2);
x(:,1) = trajTotal(:,1);
u = zeros(3,NSim);


for k = 1:NSim2-1
    % prepare variables:
    ua = -K*(x(:,k)-trajTotal(:,k));
    
    
    % compute state derivative:
    dpdt = x(4:6,k);
    dvdt = ua - Cd * x(4:6,k);

    dxdt = [dpdt; dvdt];
    
    % integrate state
    x(:,k+1) = x(:,k) + dt2 * dxdt;
end
x(:,end) = trajTotal(:,end);

% plot
figure
view(3);
hold on
plot3(trajTotal(1,:), trajTotal(2,:), trajTotal(3,:), 'LineWidth', 2);  % Start → Goal1 → Goal2
plot3(x(1,:), x(2,:), x(3,:), 'k', 'LineWidth', 2);


% Plot dos obstaculos
plot3(PStart(1), PStart(2), PStart(3), 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g')
plot3(PGoal1(1), PGoal1(2), PGoal1(3), 'mo', 'MarkerSize', 8, 'MarkerFaceColor', 'm')
plot3(PGoal2(1), PGoal2(2), PGoal2(3), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r')

axis([PFly(1,:) PFly(2,:) PFly(3,:)])

patch(obsgen(PNoFly1))
patch(obsgen(PNoFly2))
  

legend('Planeamento de Trajetória RRT*', 'Controlador desenvolvido');
title('Algoritmo de planeamento com controlador desenvolvido')



%% function

function plotPath(path, color)
    for i = 1:size(path,1)-1
        line([path(i,1), path(i+1,1)], ...
             [path(i,2), path(i+1,2)], ...
             [path(i,3), path(i+1,3)], ...
             'Color', color, 'LineWidth', 2);
    end
end
