
theta = -60:60;
for i=1:length(theta)
[CB_x(i), CB_y(i)] = get_centroid(theta(i), p.L, p.H, p.h0);
end

figure
plot(theta,CB_x-2.5,theta,CB_y)