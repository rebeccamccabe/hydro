clear;clc;close all

% inputs
theta0 = deg2rad(89);
thetadot0 = 0;
t_final = 10;
dt = 0.01;
p = struct('m', 5, 'g', 9.8, 'rho', 1000, 'Cd', .01, 'L', 5, 'H', 3, 'h0', 1, 'W', 1);
%%
% calculation
p.I = 1/12 * p.m * (p.L^2 + p.H^2);
func = @(t,x) rectDynamics(t,x,p);
sol = ode45(func,[0 t_final],[theta0 thetadot0]);

% post processing
time = 0:dt:t_final;
state = deval(sol,time);
theta = state(1,:);
thetadot = state(2,:);

% calculate vertices locations over time
ul_0 = [-p.L/2  p.H/2]'; % upper left
ur_0 = [ p.L/2  p.H/2]'; % upper right
ll_0 = [-p.L/2 -p.H/2]'; % lower left
lr_0 = [ p.L/2 -p.H/2]'; % lower right
points_0 = [ul_0 ur_0 lr_0 ll_0];

points = zeros([2 4 length(time)]);
CB = zeros([2 length(time)]);
CG = zeros([2 length(time)]);
CG_0 = [0; p.H/2-p.h0];
for i=1:length(time)
    R = [cos(theta(i)) -sin(theta(i));
        sin(theta(i)) cos(theta(i))]; % rotation matrix
    CG(:,i) = CG_0 + [0;1] .* (eye(2) - R) * CG_0; % vertical translation to maintain constant submerged area
    offset = repmat(CG(:,i), 1, 4);
    points(:,:,i) = R * points_0 + offset;
    [CB_x, CB_y] = get_centroid(rad2deg(-theta(i)), p.L, p.H, p.h0);
    CB_pre_rotate = [-CB_x; CB_y] + points_0(:,3);
    CB(:,i) = R * CB_pre_rotate + CG(:,i);
end   
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
%%
% animation plot
figure
hold on
axis equal
xlim([-p.L p.L])
ylim([-p.H 2*p.H])
line([-p.L p.L],[0 0],'HandleVisibility','off') % water line
plot(100,100,'r.','MarkerSize',20) % CG (just for legend)
plot(100,100,'k.','MarkerSize',20); % CB (just for legend)
legend('Center of Gravity','Center of Buoyancy')

for i=1:length(time)
    square = points(:,:,i);
    square_x = [square(1,:) square(1,1)];
    square_y = [square(2,:) square(2,1)];
    
    corner_x = square(1,3);
    corner_y = square(2,3);
    
    CB_x = CB(1,i);
    CB_y = CB(2,i);
    
    CG_x = CG(1,i);
    CG_y = CG(2,i);
    
    j=1:i;
    CB_x_tail = CB(1,j);
    CB_y_tail = CB(2,j);
    CG_x_tail = CG(1,j);
    CG_y_tail = CG(2,j);
    corner_x_tail = points(1,3,j);
    corner_y_tail = points(2,3,j);
    
    h1 = plot(CB_x, CB_y, 'k.',...
              CG_x, CG_y, 'r.',...
              square_x,square_y,'b-',...
              corner_x,corner_y,'b.',...
              CB_x_tail,CB_y_tail,'k',...
              CG_x_tail,CG_y_tail,'r',...
              corner_x_tail(:),corner_y_tail(:),'b',...
              'MarkerSize',20,...
              'HandleVisibility','off');
          
    %pause(3*dt)
    movieVector(i) = getframe;
    if i==length(time)
        break
    else
        delete(h1);
    end
end


makeVideo('HydroRectangle',dt,movieVector);

% dynamics function
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

function lever_arm = get_lever_arm(theta,p)
    CG = [0; p.H/2-p.h0];
    lr_0 = [p.L/2 -p.H/2]';
    R = [cos(theta) -sin(theta);
    	 sin(theta)  cos(theta)]; % rotation matrix
    [CB_x, CB_y] = get_centroid(rad2deg(-theta), p.L, p.H, p.h0);
    CB_pre_rotate = [-CB_x; CB_y] + lr_0;
    CB_rotated = R * CB_pre_rotate + CG;
    lever_arm = CG(1) - CB_rotated(1);
end