
function [C_x, C_y] = get_centroid(theta, L, H, h0)
    if theta<0
        negative = true;
        theta = abs(theta);
    else
        negative = false;
    end
    if h0 < H/2
        shape_case = 1;
        theta_c1 = atand(2 * h0 / L);
        theta_c2 = atand(H^2 / (2 * L * h0));
    elseif h0 > H/2
        shape_case = 2;
    else
        error('h0 = H/2 not coded yet')
    end

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
    
    if negative
        C_x = L - C_x;
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