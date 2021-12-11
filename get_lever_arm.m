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