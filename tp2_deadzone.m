close all;
clear; 

% Parameters
m = 2;  % kg
M = 8;  % kg
l = 0.5;  % meter
a = 1 / (m + M);

% Control gains
Kp = 1000;
Kd = 200;


A_mf = [0 1; -Kp -Kd];
B_mf = [0; 1];
btp = [0.0004 0.0016];

% Initial conditions for the pendulum
theta0 = pi/2;  % initial angle 
omega0 = 0;  % initial angular velocity

% Initial conditions for the adaptive controller
r = 100;
phi0 = zeros(r,1);

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
[t, x] = ode45(@(t, x) combined_ode(t, x, Kp, Kd, a, l, m, r, btp, A, B, C, D), tspan, [theta0; omega0; filter_state0; phi0]);


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
plot(x(:,1),x(:,2));
xlabel('\theta (rad)');
ylabel('\omega (rad/s)');
title('Trajetória no Espaço de Estados');



% Define the combined system of ODEs (filter + pendulum)
function dxdt = combined_ode(t, x, Kp, Kd, a, l, m, r, btp, A, B, C, D)
    % State variables for the filter
    filter_state = x(3:4);
    
    % State variables for the pendulum
    theta = x(1);
    omega = x(2);

    % State variable for TS Fuzzy model
    phi = x(5:end);
    
    % Reference signal
    amp = pi/2;
    freq = 1;
    theta_r = amp*sin(freq*t);
    %theta_r = pi/2 * (t >= 0);  % Unit step at t = 0

    % Nonlinear Gravity
    nu = 10;
    g0 = 9.81;
    g = g0 - nu*(cos(theta) + 0.5*sin(theta));

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

    [g_est,dphi] = approx_grav(a,m,l,r, btp, theta, omega, e_theta, e_omega, phi);
    
    den = (4*l/3) - a * m * l * cos(theta)^2;
    u = (g_est * sin(theta) - ((a * m * l * omega^2 * sin(2*theta)) / 2) - den * (-Kp * e_theta - Kd * e_omega + ddot_theta_d)) / (a * cos(theta));
    
    dtheta_dt = omega;
    domega_dt = (1 / den) * (g * sin(theta) - ((a * m * l * omega^2 * sin(2*theta)) / 2) - a * cos(theta) * u);
    
    dxdt = [dtheta_dt; domega_dt; dfilter_state_dt; dphi];
end

function [g_est,dphi] = approx_grav(a,m,l, r, btp, theta, omega, e_theta, e_omega, phi)

    g_est = 0;
    dphi = zeros(r,1);

    Gamma = 0.001;

    % Dead-zone
    delta = 0.001;

    e = [e_theta; e_omega];
    
    zeta = pertinencias(theta,r);
    
    g_est = zeta'*phi;
    
    if norm(btp*e,2) > delta
        dphi = (1/(4*l/3 - a*m*l*(cos(theta)^2)))*Gamma*zeta*(btp*e);
    else 
        dphi = zeros(r,1);
    end
end


function zeta = pertinencias(theta,r)

    lim_inf = -2;
    lim_sup = 2;

    zeta = zeros(r,1);
    if theta < lim_inf
        zeta(1) = 1;
    elseif theta > lim_sup
        zeta(end) = 1;
    else
        Theta = linspace(lim_inf,lim_sup,r); 
        zeta = interp1(Theta,eye(r),theta)';
    end
end
