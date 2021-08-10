L = 5;
H = 3;
h0 = 1;

% origin: bottom right corner of rectangle
CG_x = -L/2;
CG_y = H/2;

if h0 < H/2
    shape_case = 1;
    theta_c = atand(2 * h0 / L);
    theta_c2 = atand(H^2 / (2 * L * h0));
elseif h0 > H/2
    shape_case = 2;
else
    error('h0 = H/2 not coded yet')
end

thetas = 0:5:90;
for i=1:length(thetas)
    theta = thetas(i);
    
    if theta < theta_c_1
        % A
        Area = .5 * L *
    elseif theta == theta_c_1
        % B
        switch shape_case
            case 1
                % formula
            case 2
                % formula
        end
    elseif theta < theta_c_2
        % C
        switch shape_case
            case 1
                % formula
            case 2
                % formula
        end
    elseif theta == theta_c_2
        % D
        switch shape_case
            case 1
                % formula
            case 2
                % formula
        end
    elseif theta <= 90
        % E
        switch shape_case
            case 1
                % formula
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