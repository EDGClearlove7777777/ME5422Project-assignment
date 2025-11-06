%% ME5422 Project: Digital Controller Design for Inverted Pendulum
%
%  Author:      [Voleo Lee]
%  Date:        [05/11/2025]
%  Description: This script designs, simulates, and compares two fully
%               discrete controllers (LQR and DSMC) for stabilizing an
%               inverted pendulum.
%
%--------------------------------------------------------------------------

%% 1. Initialization and System Parameters
clear; clc; close all;

% System Physical Parameters
m0 = 0.1;   % Mass of the pendulum (Kg)
m  = 1.1;   % Mass of the cart (Kg)
L  = 0.5;   % Length to pendulum center of mass (m)
g  = 9.81;  % Acceleration due to gravity (m/s^2)

% Physical Constraints
f_max = 10;     % Actuator saturation limit (N)
x_max = 0.25;   % Cart position sensor range (m)
theta_max = 0.21; % Pendulum angle sensor range (rad)

% Simulation Parameters
Ts = 0.02;      % Sampling time (s)
T_sim = 10;     % Total simulation time (s)
N = T_sim / Ts;
time = 0:Ts:T_sim; % Time vector

% Reference Signals
x_ref_lqr = 0.15;   
x_ref_smc = -0.15;  

% Full state reference vectors (velocities are zero at steady state)
x_ref_state_lqr = [x_ref_lqr; 0; 0; 0]; 
x_ref_state_smc = [x_ref_smc; 0; 0; 0]; 

%% 2. System Modeling (Continuous-Time)

M = m + m0;
A = [0, 1, 0, 0;
     0, 0, -m0*g/m, 0;
     0, 0, 0, 1;
     0, 0, M*g/(m*L), 0];

B = [0; 1/m; 0; -1/(m*L)];

C = [1, 0, 0, 0;    % Measurable outputs: x and theta
     0, 0, 1, 0];
 
D = [0; 0];

sys_c = ss(A, B, C, D);

%% 3. System Discretization
sys_d = c2d(sys_c, Ts, 'zoh');
Ad = sys_d.A;
Bd = sys_d.B;
Cd = sys_d.C;
Dd = sys_d.D;
s_dim = size(Ad,1);

%% 4. Controller Design #1: Discrete LQR
fprintf('Designing LQR controller...\n');

% LQR Controller Tuning (Q weights state, R weights input)
Q = diag([10, 1, 100, 1]);  % Prioritize angle (100) and position (10)
R = 0.1;                    % Allow moderate control effort
[K_lqr_d, ~, ~] = dlqr(Ad, Bd, Q, R);

% State Observer Design
% Place observer poles to be 2x faster than controller poles
poles_controller_d = eig(Ad - Bd*K_lqr_d);
poles_observer_d = 0.5 * poles_controller_d; 
L_obs_d = place(Ad', Cd', poles_observer_d)'; % Use Ad' and Cd' for place()

% LQR Pre-compensator gain (for non-zero reference tracking)
M_ss = [Ad - eye(s_dim), Bd; 
        Cd(1,:),          Dd(1,1)];
       
N_xu_d = M_ss \ [zeros(s_dim,1); 1];
Nxd = N_xu_d(1:s_dim);
Nud = N_xu_d(s_dim+1);

Nbar_lqr_d = Nud + K_lqr_d * Nxd;

%% 5. Controller Design #2: Discrete Sliding Mode Control (DSMC)
fprintf('Designing DSMC controller...\n');

% 1. Define sliding surface using LQR gain (optimal surface)
C_s = K_lqr_d; 

% 2. Pre-calculate terms for the discrete control law
%    u(k) = (C_s*Bd)^-1 * [ C_s*(I - Ad)*e(k) - q*sat(s(k)) ]
inv_CsBd = 1 / (C_s * Bd);
Cs_I_minus_Ad = C_s * (eye(s_dim) - Ad);

% 3. Set SMC steady-state target vector
x_ss_smc = x_ref_state_smc;

% 4. DSMC Tuning Parameters
q_smc = 0.01;   % Reaching law gain (q)
Phi_smc = 0.05; % Boundary layer thickness (Phi) for chattering reduction

%% 6. Run Discrete-Time Simulation
fprintf('Running discrete simulation loop...\n');

% Initial Conditions
x0 = [0; 0; 0.1; 0];     % Plant: [x, x_dot, theta, theta_dot]' (starts tilted)
x_hat0 = [0; 0; 0; 0];   % Observer: (starts at equilibrium)

% Pre-allocate arrays for LQR
x_lqr = zeros(s_dim, N+1);
x_hat_lqr = zeros(s_dim, N+1);
u_lqr = zeros(1, N);
x_lqr(:,1) = x0;
x_hat_lqr(:,1) = x_hat0;

% Pre-allocate arrays for SMC
x_smc = zeros(s_dim, N+1);
x_hat_smc = zeros(s_dim, N+1);
u_smc = zeros(1, N);
x_smc(:,1) = x0;
x_hat_smc(:,1) = x_hat0;

% Run Discrete Simulation Loop
for k = 1:N
    % --- LQR Controller ---
    y_lqr = Cd * x_lqr(:,k);    % Measurement
    
    % Control law
    u_raw_lqr = -K_lqr_d * x_hat_lqr(:,k) + Nbar_lqr_d * x_ref_lqr;
    u_lqr(k) = max(-f_max, min(f_max, u_raw_lqr)); % Actuator saturation
    
    % Observer update
    x_hat_lqr(:,k+1) = Ad * x_hat_lqr(:,k) + Bd * u_lqr(k) + L_obs_d * (y_lqr - Cd * x_hat_lqr(:,k));
    % Plant update (simulated)
    x_lqr(:,k+1) = Ad * x_lqr(:,k) + Bd * u_lqr(k);

    
    % --- SMC Controller ---
    y_smc = Cd * x_smc(:,k);    % Measurement
    
    x_hat = x_hat_smc(:,k);
    e_hat = x_hat - x_ss_smc;   % Error e(k) = x_hat(k) - x_ref
    s = C_s * e_hat;            % Sliding surface s(k)
    
    % Saturation function (for boundary layer)
    if abs(s) > Phi_smc
        sat_s = sign(s);
    else
        sat_s = s / Phi_smc;
    end
    
    % Control law
    u_raw_smc = inv_CsBd * (Cs_I_minus_Ad * e_hat - q_smc * sat_s);
    u_smc(k) = max(-f_max, min(f_max, u_raw_smc)); % Actuator saturation
    
    % Observer update
    x_hat_smc(:,k+1) = Ad * x_hat_smc(:,k) + Bd * u_smc(k) + L_obs_d * (y_smc - Cd * x_hat_smc(:,k));
    % Plant update (simulated)
    x_smc(:,k+1) = Ad * x_smc(:,k) + Bd * u_smc(k);
end
fprintf('Simulation complete.\n');

%% 7. Results and Visualization
fprintf('Plotting results...\n');

% Plotting Colors
color_lqr = [0.39, 0.58, 0.93]; % Light Blue
color_smc = [1, 0.65, 0];       % Orange

% Plot: LQR + Observer Performance
figure('Name', 'Discrete LQR + Observer Performance', 'Position', [50, 50, 1000, 700]);
sgtitle('Discrete LQR + Observer Performance (Target = 0.15)', 'FontSize', 14, 'FontWeight', 'bold');
subplot(3,1,1);
plot(time, x_lqr(1,:), 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
plot([0, T_sim], [x_ref_lqr, x_ref_lqr], 'k:', 'LineWidth', 1); 
plot([0, T_sim], [x_max, x_max], 'g-.', 'LineWidth', 1);
plot([0, T_sim], [-x_max, -x_max], 'g-.', 'LineWidth', 1);
title('Cart Position (x) - LQR', 'FontSize', 12);
xlabel('Time (s)'); ylabel('Position (m)');
legend('Actual Position', 'Reference (0.15)', 'Sensor Limit', 'Location', 'southeast');
grid on; hold off;
subplot(3,1,2);
plot(time, x_lqr(3,:), 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
plot([0, T_sim], [theta_max, theta_max], 'g-.', 'LineWidth', 1);
plot([0, T_sim], [-theta_max, -theta_max], 'g-.', 'LineWidth', 1);
title('Pendulum Angle (\theta) - LQR', 'FontSize', 12);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('Actual Angle', 'Sensor Limit', 'Location', 'southeast');
grid on; hold off;
subplot(3,1,3);
stairs(time(1:end-1), u_lqr, 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
plot([0, T_sim], [f_max, f_max], 'k:', 'LineWidth', 1);
plot([0, T_sim], [-f_max, -f_max], 'k:', 'LineWidth', 1);
title('Control Input (Force) - LQR (Digital/Stepped)', 'FontSize', 12);
xlabel('Time (s)'); ylabel('Force (N)');
legend('Control Force', 'Saturation Limit', 'Location', 'southeast');
grid on; hold off;

% Plot: Discrete SMC Controller Performance
figure('Name', 'Discrete SMC Controller Performance', 'Position', [150, 150, 1000, 700]);
sgtitle('Fully Discrete Sliding Mode Control (DSMC) (Target = -0.15)', 'FontSize', 14, 'FontWeight', 'bold');
subplot(3,1,1);
plot(time, x_smc(1,:), 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
hold on;
plot([0, T_sim], [x_ref_smc, x_ref_smc], 'k:', 'LineWidth', 1); 
plot([0, T_sim], [x_max, x_max], 'g-.', 'LineWidth', 1);
plot([0, T_sim], [-x_max, -x_max], 'g-.', 'LineWidth', 1);
title('Cart Position (x) - DSMC', 'FontSize', 12);
xlabel('Time (s)'); ylabel('Position (m)');
legend('Actual Position', 'Reference (-0.15)', 'Sensor Limit', 'Location', 'southeast');
grid on; hold off;
subplot(3,1,2);
plot(time, x_smc(3,:), 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
hold on;
plot([0, T_sim], [theta_max, theta_max], 'g-.', 'LineWidth', 1);
plot([0, T_sim], [-theta_max, -theta_max], 'g-.', 'LineWidth', 1);
title('Pendulum Angle (\theta) - DSMC', 'FontSize', 12);
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('Actual Angle', 'Sensor Limit', 'Location', 'southeast');
grid on; hold off;
subplot(3,1,3);
stairs(time(1:end-1), u_smc, 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
hold on;
plot([0, T_sim], [f_max, f_max], 'k:', 'LineWidth', 1);
plot([0, T_sim], [-f_max, -f_max], 'k:', 'LineWidth', 1);
title('Control Input (Force) - DSMC (Digital/Stepped)', 'FontSize', 12);
xlabel('Time (s)'); ylabel('Force (N)');
legend('Control Force', 'Saturation Limit', 'Location', 'southeast');
grid on; hold off;


% Plot: Combined Comparison
figure('Name', 'Controller Performance Comparison (FULLY DIGITAL)', 'Position', [300, 100, 1200, 800]);
sgtitle('Direct Comparison: Digital LQR (0.15) vs. Discrete SMC (-0.15)', 'FontSize', 14, 'FontWeight', 'bold');
subplot(3,1,1);
plot(time, x_lqr(1,:), 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
plot(time, x_smc(1,:), 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
plot([0, T_sim], [x_ref_lqr, x_ref_lqr], ':', 'Color', color_lqr, 'LineWidth', 2);
plot([0, T_sim], [x_ref_smc, x_ref_smc], ':', 'Color', color_smc, 'LineWidth', 2);
plot([0, T_sim], [x_max, x_max], 'm-.', 'LineWidth', 1);
plot([0, T_sim], [-x_max, -x_max], 'm-.', 'LineWidth', 1);
title('Cart Position (x) Comparison');
xlabel('Time (s)'); ylabel('Position (m)');
legend('LQR (Light Blue)', 'DSMC (Orange)', 'LQR Ref (0.15)', 'SMC Ref (-0.15)', 'Sensor Limit', 'Location', 'best'); 
grid on; hold off;
subplot(3,1,2);
plot(time, x_lqr(3,:), 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
plot(time, x_smc(3,:), 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
plot([0, T_sim], [theta_max, theta_max], 'm-.', 'LineWidth', 1);
plot([0, T_sim], [-theta_max, -theta_max], 'm-.', 'LineWidth', 1);
title('Pendulum Angle (\theta) Comparison');
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('LQR (Light Blue)', 'DSMC (Orange)', 'Sensor Limit', 'Location', 'southeast'); 
grid on; hold off;
subplot(3,1,3);
stairs(time(1:end-1), u_lqr, 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
stairs(time(1:end-1), u_smc, 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
plot([0, T_sim], [f_max, f_max], 'k:', 'LineWidth', 1);
plot([0, T_sim], [-f_max, -f_max], 'k:', 'LineWidth', 1);
title('Control Input  Comparison ');
xlabel('Time (s)'); ylabel('Force (N)');
legend('LQR (Light Blue)', 'DSMC (Orange)', 'Saturation Limit', 'Location', 'southeast'); 
grid on; hold off;


% Plot: Separate Cart Position Comparison
figure('Name', 'Cart Position Comparison (Separate)', 'Position', [400, 200, 800, 500]);
plot(time, x_lqr(1,:), 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
plot(time, x_smc(1,:), 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
plot([0, T_sim], [x_ref_lqr, x_ref_lqr], ':', 'Color', color_lqr, 'LineWidth', 2);
plot([0, T_sim], [x_ref_smc, x_ref_smc], ':', 'Color', color_smc, 'LineWidth', 2);
plot([0, T_sim], [x_max, x_max], 'm-.', 'LineWidth', 1);
plot([0, T_sim], [-x_max, -x_max], 'm-.', 'LineWidth', 1);
title('Cart Position (x) Comparison: LQR (0.15) vs. DSMC (-0.15)');
xlabel('Time (s)'); ylabel('Position (m)');
legend('LQR (Light Blue)', 'DSMC (Orange)', 'LQR Ref (0.15)', 'SMC Ref (-0.15)', 'Sensor Limit', 'Location', 'best'); 
grid on; hold off;

% Plot: Separate Pendulum Angle Comparison
figure('Name', 'Pendulum Angle Comparison (Separate)', 'Position', [500, 250, 800, 500]);
plot(time, x_lqr(3,:), 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold on;
plot(time, x_smc(3,:), 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
plot([0, T_sim], [theta_max, theta_max], 'm-.', 'LineWidth', 1);
plot([0, T_sim], [-theta_max, -theta_max], 'm-.', 'LineWidth', 1);
title('Pendulum Angle (\theta) Comparison: LQR (0.15) vs. DSMC (-0.15)');
xlabel('Time (s)'); ylabel('Angle (rad)');
legend('LQR (Light Blue)', 'DSMC (Orange)', 'Sensor Limit', 'Location', 'southeast'); 
grid on; hold off;


% Plot: Control Input Comparison (with Inset)
figure('Name', 'Control Input Comparison (with Inset)', 'Position', [600, 300, 800, 500]);

% 1. Plot main axes
ax1 = gca; 
stairs(ax1, time(1:end-1), u_lqr, 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2); 
hold(ax1, 'on');
stairs(ax1, time(1:end-1), u_smc, 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2); 
plot(ax1, [0, T_sim], [f_max, f_max], 'k:', 'LineWidth', 1);
plot(ax1, [0, T_sim], [-f_max, -f_max], 'k:', 'LineWidth', 1);
title(ax1, 'Control Input  Comparison ');
xlabel(ax1, 'Time (s)'); ylabel(ax1, 'Force (N)');
legend(ax1, 'LQR (Light Blue)', 'DSMC (Orange)', 'Saturation Limit', 'Location', 'southeast'); 
grid(ax1, 'on');

% 2. Draw rectangle for zoom area
rectangle(ax1, 'Position', [0, -f_max, 1, 2*f_max], 'EdgeColor', 'r', 'LineStyle', '--', 'LineWidth', 1);

% 3. Create inset axes
%    [left, bottom, width, height] (normalized coords)
ax2 = axes('Position', [0.45, 0.4, 0.4, 0.4]); 
stairs(ax2, time(1:end-1), u_lqr, 'LineStyle', '-', 'Color', color_lqr, 'LineWidth', 2);
hold(ax2, 'on');
stairs(ax2, time(1:end-1), u_smc, 'LineStyle', '-', 'Color', color_smc, 'LineWidth', 2);
plot(ax2, [0, T_sim], [f_max, f_max], 'k:', 'LineWidth', 1);
plot(ax2, [0, T_sim], [-f_max, -f_max], 'k:', 'LineWidth', 1);

% 4. Set inset limits
xlim(ax2, [0, 1]);
ylim(ax2, [-f_max, f_max]); % Match Y-axis

% 5. Format inset
grid(ax2, 'on');
box(ax2, 'on');
set(ax2, 'XTick', 0:0.2:1); 

hold(ax1, 'off'); 

%% 8. Animation
fprintf('Generating animations...\n');
animate_pendulum(time, x_lqr, L, 'Digital LQR (Target 0.15)');
animate_pendulum(time, x_smc, L, 'Discrete SMC (Target -0.15)');


%% 9. Animation Function
function animate_pendulum(t, x_data, L, controller_name)
    % x_data = [x; x_dot; theta; theta_dot]
    
    x_pos = x_data(1,:);
    theta_pos = x_data(3,:);
    
    figure('Name', ['Animation: ' controller_name], 'Position', [randi(800), randi(400), 800, 400]);
    ax = gca;
    ax.NextPlot = 'replaceChildren';
    
    cart_width = 0.2;
    cart_height = 0.1;

    frame_skip = max(1, floor(length(t) / 200)); % Aim for ~200 frames

    for k = 1:frame_skip:length(t)
        if isnan(x_pos(k)) || isnan(theta_pos(k)) 
            break;
        end
        cla; 
        
        % Draw Cart
        cart_x = x_pos(k);
        patch(ax, [cart_x-cart_width/2, cart_x+cart_width/2, cart_x+cart_width/2, cart_x-cart_width/2], ...
                  [0, 0, cart_height, cart_height], 'blue');
        
        % Draw Pendulum
        pendulum_x = cart_x + L * sin(theta_pos(k));
        pendulum_y = cart_height + L * cos(theta_pos(k));
        line(ax, [cart_x, pendulum_x], [cart_height, pendulum_y], 'Color', 'red', 'LineWidth', 3);
        
        % Draw Pendulum Mass
        viscircles(ax, [pendulum_x, pendulum_y], 0.05, 'Color', 'black', 'LineWidth', 2);
        line(ax, [-1, 1], [0, 0], 'Color', 'black'); % Ground
        
        % Widen axis limits to see both targets
        axis(ax, [-0.6, 0.6, -0.1, 1.2*L]); 
        daspect(ax, [1,1,1]);
        title(ax, sprintf('%s\nTime: %.2f s', controller_name, t(k)));
        
        drawnow;
    end
end
