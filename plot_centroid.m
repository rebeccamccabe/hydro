CG = [0; p.H/2-p.h0];
theta = -60:60;

for i=1:length(theta)
[CB_x(i), CB_y(i)] = get_centroid(theta(i), p.L, p.H, p.h0);

R = [cosd(theta(i)) -sind(theta(i));
        sind(theta(i)) cosd(theta(i))];
translation(:,i) = (eye(2) - R) * CG;
end

figure
plot(theta,CB_x-2.5,theta,CB_y)
figure
plot(theta,translation(1,:),theta,translation(2,:))