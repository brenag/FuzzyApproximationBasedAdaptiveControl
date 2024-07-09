close all;
clc; 

% Parameters
m = 2;  % kg
M = 8;  % kg
l = 0.5;  % meter
a = 1 / (m + M);

% Control gains
Kp = 1000;
Kd = 200;

% Initial conditions for the pendulum
theta0 = pi/2;  % initial angle (45 degrees)
omega0 = 0;  % initial angular velocity

% Time span
tspan = [0 30];

% State-space model for the trajectory generation filter
A = [0 1; -100 -20];
B = [0; 100];
C = [1 0; 0 1; -100 -20];
D = [0; 0; 100];

% Initial conditions for the filter
filter_state0 = [0; 0];

% Solve the system using ode45
[t, x] = ode45(@(t, x) combined_ode(t, x, Kp, Kd, a, l, m, A, B, C, D), tspan, [ theta0; omega0; filter_state0]);


% Plot results
figure;
subplot(2, 1, 1);
plot(t, x(:, 1), t, x(:,3));
xlabel('Tempo (s)');
ylabel('\theta (rad)');
title('Posição angular \theta no Tempo');
legend('\theta', '\theta_d')

subplot(2, 1, 2);
plot(t, x(:, 2),t,x(:,4));
xlabel('Tempo (s)');
ylabel('\omega (rad/s)');
title('Velocidade Angular \omega no Tempo');
legend('\omega', '\omega_d')

figure;
subplot(2, 1, 1);
plot(t, x(:, 1) - x(:,3));
xlabel('Tempo (s)');
ylabel('e_{\theta} (rad)');
title('Erro de posição angular no Tempo');

subplot(2, 1, 2);
plot(t, x(:, 2) - x(:,4));
xlabel('Tempo (s)');
ylabel('e_{\omega} (rad/s)');
title('Erro de velocidade angular no Tempo');

% Define the combined system of ODEs (filter + pendulum)
function dxdt = combined_ode(t, x, Kp, Kd, a, l, m, A, B, C, D)
    % State variables for the filter
    filter_state = x(3:4);
    
    % State variables for the pendulum
    theta = x(1);
    omega = x(2);
    
    % Reference signal
    theta_r = pi/2 *sin(t);
    %theta_r = pi/2 * (t >= 0);  % Unit step at t = 0

    % Trajectory generation filter (state-space model)
    dfilter_state_dt = A * filter_state + B * theta_r;
    filtered_output = C * filter_state + D * theta_r;
    
    % Filtered references
    theta_d = filtered_output(1);
    dot_theta_d = filtered_output(2);
    ddot_theta_d = filtered_output(3);
    
    % Pendulum control
    e_theta = theta - theta_d;
    e_omega = omega - dot_theta_d;
    
    nu = 10;
    g0 = 9.81;
    g = g0 - nu*(cos(theta) + 0.5*sin(theta));
    
    den = (4*l/3) - a * m * l * cos(theta)^2;
    u = (g * sin(theta) - (a * m * l * omega^2 * sin(2*theta)) / 2 - den * (-Kp * e_theta - Kd * e_omega + ddot_theta_d)) / (a * cos(theta));
    
    dtheta_dt = omega;
    domega_dt = (1 / den) * (g * sin(theta) - ((a * m * l * omega^2 * sin(2*theta)) / 2) - a * cos(theta) * u);
    
    dxdt = [dtheta_dt; domega_dt; dfilter_state_dt];
end
