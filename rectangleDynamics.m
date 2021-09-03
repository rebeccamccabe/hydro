clear;clc;close all

% inputs
theta0 = deg2rad(45);
thetadot0 = 0;
t_final = 10;
dt = 0.01;
p = struct('m', 5, 'g', 9.8, 'rho', 1000, 'Cd', .01, 'L', 5, 'H', 3, 'h0', 1, 'W', 1);

% calculation
p.I = 1/12 * p.m * (p.L^2 + p.H^2);
func = @(t,x)rectDynamics(t,x,p);
sol = ode45(func,[0 t_final],[theta0 thetadot0]);

% post processing
time = 0:dt:t_final;
state = deval(sol,time);
theta = state(1,:);
thetadot = state(2,:);

% time series plot
figure
plot(time,theta,time,thetadot)
legend('theta vs time', 'thetadot vs time')

% calculate vertices locations over time
ul_0 = [-p.L/2  p.H/2]'; % upper left
ur_0 = [ p.L/2  p.H/2]'; % upper right
ll_0 = [-p.L/2 -p.H/2]'; % lower left
lr_0 = [ p.L/2 -p.H/2]'; % lower right
points_0 = [ul_0 ur_0 lr_0 ll_0];

points = zeros([2 4 length(time)]);
CB_world = zeros([2 length(time)]);
CG = [0; p.H/2-p.h0];
for i=1:length(time)
    R = [cos(theta(i)) -sin(theta(i));
        sin(theta(i)) cos(theta(i))]; % rotation matrix
    offset = repmat(CG, 1, 4);
    points(:,:,i) = R * points_0 + offset;
    
    [CB_x, CB_y] = get_centroid(theta(i), p.L, p.H, p.h0);
    CB_world(:,i) = [1 -1; -1 1] .* R * [CB_x; CB_y] .* [-1; 1] + points(:,3,i);
end   

% animation plot
figure
hold on
axis equal
xlim([-2*p.L 2*p.L])
ylim([-2*p.H 2*p.H])
line([-2*p.L 2*p.L],[0 0]) % water line
plot(CG(1),CG(2),'ko') % center of gravity

for i=1:length(time)
    square = points(:,:,i);
    square_x = [square(1,:) square(1,1)];
    square_y = [square(2,:) square(2,1)];
    CBx = CB_world(1,i);
    CBy = CB_world(2,i);
    
    h1 = plot(square_x,square_y,'b-');%,CBx, CBy,'ko');
    %patch(square(1,:),square(2,:),'b');
    %h1 = plot(polyshape(square'),'FaceColor','b');
    
    pause(dt)
    movieVector(i) = getframe;
    if i==length(time)
        break
    else
        delete(h1);
    end
end

%%
myWriter = VideoWriter('HydroRectangle','MPEG-4');
myWriter.FrameRate = 1/dt;
open(myWriter);
writeVideo(myWriter,movieVector);
close(myWriter);

% dynamics function
function dxdt = rectDynamics(~,x,p)

theta = x(1);
thetadot = x(2);

lever_arm = sin(theta);
tau_b = -p.m * p.g * lever_arm; % buoyancy torque
tau_d = -1/2 * p.rho * p.Cd * p.W * p.h0^4 * thetadot * abs(thetadot); % damping torque

dxdt = [0;0];
dxdt(1) = thetadot;
dxdt(2) = 1/p.I * (tau_b + tau_d);

end