clear;clc;close all

% inputs
theta0 = deg2rad(89);
thetadot0 = 0;
t_final = 10;
dt = 0.01;
p = parameters();

%% calculation
[state,time] = ode45(@rectDynamics,t_final,[theta0 thetadot0],p,dt);

%% post processing
theta = state(1,:);
thetadot = state(2,:);

size = length(time);
[points, CG, CB] = get_points(p,size,theta);

%%
% time series plot
y = CG(2,:);
figure
plot(time,theta,time,thetadot,time,y)
legend('theta vs time', 'thetadot vs time', 'CG height vs time')
%%
figure
yddot = [0 diff(diff(y)) 0];
deltaFb = p.m*yddot;
percentDeltaFb = 100 * deltaFb / (p.m*p.g);
plot(time,percentDeltaFb)

%% animation
saveMovie = true;
animate(p,time,points,CB,CG,saveMovie)

%% dynamics function
function dxdt = rectDynamics(~,x,p)

theta = x(1);
thetadot = x(2);

lever_arm = get_lever_arm(theta, p); % center of buoyancy lever arm using centroid math
%lever_arm = sin(theta); % for comparison to pendulum
tau_b = -p.m * p.g * lever_arm; % buoyancy torque
tau_d = -1/2 * p.rho * p.Cd * p.W * p.h0^4 * thetadot * abs(thetadot); % damping torque

dxdt = [0;0];
dxdt(1) = thetadot;
dxdt(2) = 1/p.I * (tau_b + tau_d);

end

