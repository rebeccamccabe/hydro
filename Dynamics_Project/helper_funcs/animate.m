function [] = animate(p,time,points,CB,CG,saveMovie,plot_title)

dt = time(2)-time(1);

figure
hold on
axis equal

x = linspace(-2*p.L, 2*p.L, 100);

xlim([x(1) x(end)])
ylim([-p.H 2*p.H])

xlabel('x')
ylabel('y')
title(plot_title)

line([x(1) x(end)],[0 0],'HandleVisibility','off') % water line
wavelength = p.g*pi/p.w^2;
plot(x,p.Hs*sin(2*pi/wavelength*x),'b','HandleVisibility','off') % wave
plot(0,-p.h0,'kx','HandleVisibility','off') % pivot point

plot(100,100,'r.','MarkerSize',20) % CG (just for legend)
plot(100,100,'k.','MarkerSize',20); % CB (just for legend)
legend('Center of Gravity','Center of Buoyancy')

for i=1:length(time)
    square = points(:,1:4,i);
    square_x = [square(1,:) square(1,1)];
    square_y = [square(2,:) square(2,1)];
    
    if size(points,2) > 4
        square2 = points(:,[1 2 6 5],i);
        square2_x = [square2(1,:) square2(1,1)];
        square2_y = [square2(2,:) square2(2,1)];
        idx = 5;
    else
        square2_x = NaN;
        square2_y = NaN;
        idx = 3;
    end
    
    corner_x = points(1,idx,i);
    corner_y = points(2,idx,i);
    
    CB_x = CB(1,i);
    CB_y = CB(2,i);
    
    CG_x = CG(1,i);
    CG_y = CG(2,i);
    
    j=1:i;
    CB_x_tail = CB(1,j);
    CB_y_tail = CB(2,j);
    CG_x_tail = CG(1,j);
    CG_y_tail = CG(2,j);
    corner_x_tail = points(1,idx,j);
    corner_y_tail = points(2,idx,j);
    
    h1 = plot(CB_x, CB_y, 'k.',...
              CG_x, CG_y, 'r.',...
              square_x,square_y,'b-',...
              square2_x,square2_y,'g-',...
              corner_x,corner_y,'g.',...
              CB_x_tail,CB_y_tail,'k',...
              CG_x_tail,CG_y_tail,'r',...
              corner_x_tail(:),corner_y_tail(:),'g',...
              'MarkerSize',20,...
              'HandleVisibility','off');
          
    if saveMovie
        movieVector(i) = getframe;
    else
        pause(3*dt)
    end
    
    if i==length(time)
        break
    else
        delete(h1);
    end
end

if saveMovie
    makeVideo('HydroRectangle',dt,movieVector);
end

end