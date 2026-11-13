%% Description
%{
Main script for analyzing a ballistic hypersonic vehicle trajectory. Calculates peak 
aerodynamic heating and deceleration using analytical expressions, then compares the 
results against numerical simulations using both simplified ballistic and full 
trajectory equations. Generates a comparison table with percent errors and plots 
velocity and time as functions of altitude.
%}
%% Givens
g = 9.81; % [m/s^2]
R_e = 6378e3; % radius of Earth [m]
rho_SL = 1.225; % [kg/m^3]
alpha = 1.378e-4; % [1/m]

R_A = 0.05; % nose radius [m]
R_B = 0.9; % base radius [m]
A = pi * R_B^2;
m = 3400; % [kg]
C_D = 0.4;
C_L = 0;
L = 0;
T = 0;

K = 2e-4/sqrt(R_A); % [kg^1/2 / m]
V0 = 4500; % [m/s]
h0 = 100000; % [m]
gamma0 = deg2rad(-70); % [rad]

dt = 0.001; % [s]
W = m*g;

%% Analytical Expressions
beta = -2*W*alpha*sin(gamma0) / (g*C_D*A); % [kg/m^3]

% Deceleration
V_a_max = V0*exp(-1/2);
V_dot_max = V0^2*alpha*sin(gamma0)/(2*exp(1)*g);
h_a_max = -1/alpha * log(-beta/rho_SL * log(V_a_max/V0));

% Heating
V_q_max = V0*exp(-1/6);
q_max = K*sqrt((-beta*log(V_q_max/V0)))*V_q_max^3;
h_q_max = -1/alpha * log(-beta/rho_SL * log(V_q_max/V0));

%% Simplified Ballistic
[t_ball, h_ball, V_ball, gamma_ball, q_dot_ball, accel_ball, i_ball] = ...
    Trajectory('ballistic', V0, gamma0, h0, dt, g, R_e, rho_SL, alpha, A, C_D, C_L, T, K, W);

% Deceleration
[V_dot_max_ball, i_V_dot_max_ball] = min(accel_ball);
V_max_ball = V_ball(i_V_dot_max_ball);
h_max_ball = h_ball(i_V_dot_max_ball);
t_max_ball = t_ball(i_V_dot_max_ball);

% Heating
[q_max_ball, i_q_max_ball] = max(q_dot_ball);
V_q_max_ball = V_ball(i_q_max_ball);
h_q_max_ball = h_ball(i_q_max_ball);
t_q_max_ball = t_ball(i_q_max_ball);

%% Full Trajectory
[t_full, h_full, V_full, gamma_full, q_dot_full, accel_full, i_full] = ...
    Trajectory('full', V0, gamma0, h0, dt, g, R_e, rho_SL, alpha, A, C_D, C_L, T, K, W);

% Deceleration
[V_dot_max_full, i_V_dot_max_full] = min(accel_full);
V_max_full = V_full(i_V_dot_max_full);
h_max_full = h_full(i_V_dot_max_full);
t_max_full = t_full(i_V_dot_max_full);

pct_V_dot_full = abs(100*(V_dot_max_full - V_dot_max)/V_dot_max);
pct_V_Vdot_full = abs(100*(V_max_full - V_a_max)/V_a_max);
pct_h_Vdot_full = abs(100*(h_max_full - h_a_max)/h_a_max);

% Heating
[q_max_full, i_q_max_full] = max(q_dot_full);
V_q_max_full = V_full(i_q_max_full);
h_q_max_full = h_full(i_q_max_full);
t_q_max_full = t_full(i_q_max_full);

pct_q_max_full = abs(100*(q_max_full - q_max)/q_max);
pct_h_q_full = abs(100*(h_q_max_full - h_q_max)/h_q_max);
pct_V_q_full = abs(100*(V_q_max_full - V_q_max)/V_q_max);

%% Results
RowNames = {'Peak Deceleration [g]'; 'V @ Peak Decel [m/s]'; 'h @ Peak Decel [km]'; ...
            'Peak Heat Flux [W/m^2]'; 'V @ Peak Heat Flux [m/s]'; 'h @ Peak Heat Flux [km]'};
 
Analytical = [V_dot_max;      V_a_max;      h_a_max/1000; ...
              q_max;      V_q_max;      h_q_max/1000];
 
Ballistic  = [V_dot_max_ball; V_max_ball;   h_max_ball/1000; ...
              q_max_ball; V_q_max_ball; h_q_max_ball/1000];
 
Full       = [V_dot_max_full; V_max_full;   h_max_full/1000; ...
              q_max_full; V_q_max_full; h_q_max_full/1000];
 
Error_Full_vs_Simplification = [pct_V_dot_full;  pct_V_Vdot_full; pct_h_Vdot_full; ...
                pct_q_max_full;  pct_V_q_full;    pct_h_q_full];
 
ResultsTable = table(Analytical, Ballistic, Full, Error_Full_vs_Simplification, ...
    'RowNames', RowNames)

%% Plot
figure('Name', 'Velocity and Time vs Altitude');

subplot(1,2,1);
plot(V_full/V0, h_full/1000);
hold on;
plot(V_ball/V0, h_ball/1000);
xlabel('V/V0');
ylabel('Altitude, h [km]');
title('Velocity vs Altitude');
legend('Full Trajectory', 'Ballistic Simplification', 'Location', 'best');
hold off;

subplot(1,2,2)
plot(t_full, h_full/1000);
hold on;
plot(t_ball, h_ball/1000);
xlabel('Time, t [s]');
ylabel('Altitude, h [km]');
title('Time vs Altitude');
legend('Full Trajectory', 'Ballistic Simplication', 'Location', 'best');
hold off;

figure(2);
plot(rad2deg(gamma_full), h_full/1000);
xlabel('gamma [deg]');
ylabel('Altitude, h [km]');
%% Functions Used
function rho = density(h, rho_SL, alpha)
    rho = rho_SL*exp(-alpha*h);
end

function D = drag(rho, V, C_D, S)
    D = 0.5*rho*V^2*C_D*S;
end

function L = lift(rho, V, C_L, S)
    L = 0.5*rho*V^2*C_L*S;
end

function qdot = heatflux(rho, V, K)
    qdot = K*sqrt(rho)*V^3;
end

