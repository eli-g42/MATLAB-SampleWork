%% Description
%{
Function that numerically integrates hypersonic vehicle trajectory equations using 
first-order explicit Euler time integration with a 1 ms time step. Supports both 
full trajectory equations and simplified ballistic assumptions, calculating velocity, 
flight path angle, altitude, aerodynamic heating, and deceleration throughout the 
vehicle's descent.
%}


function [t, h, V, gamma, q_dot, accel, i] = ...
    Trajectory(mode, V0, gamma0, h0, dt, g, R_e, rho_SL, alpha, A, C_D, C_L, T, K, W)

% mode = 'full'      : integrate full ODE system, general L and T
% mode = 'ballistic' : gamma constant, gravity term dropped from
%                      V-dot (valid when W*sin(gamma) << D), L = T = 0 assumed
i = 0;
% t_max = 300; % s, safety cutoff

%% Density, Drag, Lift, and heat flux (q-dot)
density = @(h) rho_SL * exp(-alpha*h);
drag = @(rho, V) 0.5*rho*V^2*C_D*A;
lift = @(rho,V) 0.5*rho*V^2*C_L*A;
heat_flux = @(rho, V) K*sqrt(rho)*V^3;

%% Initial Conditions
t(1) = 0;
h(1) = h0;
V(1) = V0;
gamma(1) = gamma0;

%% 
while h(i+1) > 0 % && t(i+1) < t_max
    R   = R_e + h(i+1);
    rho = density(h(i+1));
    D   = drag(rho, V(i+1));
    L   = lift(rho, V(i+1));

    switch mode
        case 'full'
            V_dot = -g*D/W - g*sin(gamma(i+1));
            gamma_dot = 1/V(i+1) * (g*(L/W) - (g - V(i+1)^2/R)*cos(gamma(i+1)));

        case 'ballistic'
            V_dot = -g*D/W;
            gamma_dot = 0;
    end

    q_dot(i+1) = heat_flux(rho, V(i+1));
    accel(i+1) = V_dot/g;

    % Euler update
    V(i+2) = V(i+1) + V_dot*dt;
    gamma(i+2) = gamma(i+1) + gamma_dot*dt;
    h(i+2) = h(i+1) + V(i+1)*sin(gamma(i+1))*dt; % note: V*sin(gamma) = h_dot
    t(i+2) = t(i+1) + dt;

    i = i+1; 
end

t = t(1:end-1);
h = h(1:end-1);
V = V(1:end-1);
gamma = gamma(1:end-1);
end
