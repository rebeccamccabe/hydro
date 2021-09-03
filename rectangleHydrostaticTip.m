clear;close all;clc

%% input parameters
L = 5;
H = 3;
h0 = 1;
CG_x_offset_ratio = 0.20; % 0 is no offset, .5 is fully offset

%% calculations

% origin: bottom right corner of rectangle
% x positive is left
% y positive is up
CG_x = L/2 - CG_x_offset_ratio*L;
CG_y = H/2;

thetas = 0:5:90;
tipping = logical(size(thetas));
angles = zeros(size(thetas));
for i=1:length(thetas)
    theta = thetas(i);
    [CB_x, CB_y] = get_centroid(theta, L, H, h0);
    [tipping(i), angles(i)] = see_if_tipping(theta, CB_x, CB_y, CG_x, CG_y);
end

figure
plot(thetas, tipping, thetas, angles)
xlabel('Degrees')
ylabel('Tipping?')
%improvePlot

%% helper functions
function [tipping, angle] = see_if_tipping(theta, CB_x, CB_y, CG_x, CG_y)
    tipping = false;
    
    delta_x = CG_x - CB_x;
    delta_y = CG_y - CB_y;
    alpha = atand(delta_x / delta_y);
    angle = alpha - theta;
    
    if angle < 0
        tipping = true;
    end
    if any([delta_x delta_y]<0)
        %tipping = NaN;
    end
end
