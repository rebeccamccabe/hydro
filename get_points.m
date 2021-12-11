function [points, CG, CB] = get_points(p,num_pts,theta,r)

    % calculate vertices locations over time
    if nargin>3
        ul_0 = [-p.L/2  p.H]'; % upper left
        ur_0 = [ p.L/2  p.H]'; % upper right
        ll_0 = [-p.L/2  0]';   % lower left
        lr_0 = [ p.L/2  0]';   % lower right
        ul2_0 = ul_0;          % upper left extension
        ur2_0 = ur_0;          % upper right extension
        points_0 = [ul_0 ur_0 lr_0 ll_0 ul2_0 ur2_0];
        CG_0 = [0; p.H/2];
    else
        ul_0 = [-p.L/2  p.H/2]'; % upper left
        ur_0 = [ p.L/2  p.H/2]'; % upper right
        ll_0 = [-p.L/2 -p.H/2]'; % lower left
        lr_0 = [ p.L/2 -p.H/2]'; % lower right
        points_0 = [ul_0 ur_0 lr_0 ll_0];
        CG_0 = [0; p.H/2-p.h0];
    end
    

    points = zeros([2 size(points_0,2) num_pts]);
    CB = zeros([2 num_pts]);
    CG = zeros([2 num_pts]);
    
    
    for i = 1:num_pts
        R = [cos(theta(i)) -sin(theta(i));
            sin(theta(i)) cos(theta(i))]; % rotation matrix
        
        if nargin>3 % r dynamics are actually simulated
            offset = [0; -p.h0];
            CG(:,i) = CG_0' + r(i) .* R(2,:) + offset';
            points_adjust = [0 0 0 0 0    0;     % x
                             0 0 0 0 r(i) r(i)]; % y
        else % r is calculated from results assuming constant buoyancy
            CG(:,i) = CG_0 + [0;1] .* (eye(2) - R) * CG_0; % vertical translation to maintain constant submerged area
            offset = CG(:,i);
            points_adjust = 0;
        end
        [CB_x, CB_y] = get_centroid(rad2deg(-theta(i)), p.L, p.H, p.h0);
        CB_pre_rotate = [-CB_x; CB_y] + points_0(:,3);
        CB(:,i) = R * CB_pre_rotate + offset;
        points(:,:,i) = R * (points_0 + points_adjust) + repmat(offset,1,size(points_0,2));
    end   
    
end
