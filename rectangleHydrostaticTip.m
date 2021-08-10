clear;close all;clc

%% input parameters
L = 5;
H = 3;
h0 = 1;
CG_x_offset_ratio = 0.20;

%% calculations

% origin: bottom right corner of rectangle
% x positive is left
% y positive is up
CG_x = L/2 - CG_x_offset_ratio*L;
CG_y = H/2;

if h0 < H/2
    shape_case = 1;
    theta_c1 = atand(2 * h0 / L);
    theta_c2 = atand(H^2 / (2 * L * h0));
elseif h0 > H/2
    shape_case = 2;
else
    error('h0 = H/2 not coded yet')
end

thetas = 0:5:90;
tipping = logical(size(thetas));
angles = zeros(size(thetas));
for i=1:length(thetas)
    theta = thetas(i);
    [CB_x, CB_y] = get_centroid(theta, theta_c1, theta_c2, L, H, h0, shape_case);
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

function [C_x, C_y] = get_centroid(theta, theta_c1, theta_c2, L, H, h0, shape_case)
    if theta < theta_c1
        % A
        h2 = h0 + L/2*tand(theta);
        h1 = h0 - L/2*tand(theta);
        [C_x, C_y] = centroid_trapezoid(L, h2, h1);
%     elseif theta == theta_c1
%         % B
%         switch shape_case
%             case 1
%                 % formula
%             case 2
%                 % formula
%         end
    elseif theta <= theta_c2
        % B, C, D
        switch shape_case
            case 1
                l1 = sqrt(2*L*h0/tand(theta));
                h1 = sqrt(2*L*h0*tand(theta));
                [C_x, C_y] = centroid_triangle( l1, h1 );
            case 2
                % formula
        end
%     elseif theta == theta_c2
%         % D
%         switch shape_case
%             case 1
%                 % formula
%             case 2
%                 % formula
%         end
    elseif theta <= 90
        % E
        switch shape_case
            case 1
                l2 = L*h0/H - H/(2*tand(theta));
                l1 = L*h0/H + H/(2*tand(theta));
                [C_x_tmp, C_y_tmp] = centroid_trapezoid(H, l2, l1);
                C_x = C_y_tmp;
                C_y = H - C_x_tmp;
            case 2
                % formula
        end
    else
        error('theta > 90 not supported')
    end
end

function [C_x, C_y] = centroid_triangle(base, height)
%                  /|
%                /  |
%              /    | height   
%            /   x  |               ^ 
%          /        |               | C_y
%         ----------Origin          |
%             base
%                 <-- C_x
%
    C_x = base/3;
    C_y = height/3;
end

function [C_x, C_y] = centroid_trapezoid(base, side_right, side_left)
%                  /|
%                /  |
%              /    |
%            /      |
% side_left |    x  | side_right    ^
%           |       |               | C_y
%           --------Origin          |
%             base
%                 <-- C_x
%
% equation source: https://www.efunda.com/math/areas/Trapezoid.cfm

    h = base;
    a = side_left;
    b = side_right;
    
    C_x = h*(2*a+b) /(3*(a+b));
    C_y = (a^2 + a*b + b^2) /(3*(a+b));
end