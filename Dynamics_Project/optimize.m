% trying to find the optimal stiffness and damping of the powertrain
% the optimization doesn't seem to work very well
% note: x = [bx bt kx kt];

opts = optimoptions('fmincon','PlotFcn',@optimplotfval,'OptimalityTolerance',1e-3);
min = [1 1 1 1]*-1e4;
max = [1 1 1 1]*1e4;
x0 = [1 1 1 1]; 
[x_opt,neg_power] = fmincon(@(x)hydro_wrap(x,false),x0,[],[],[],[],min,max,[],opts);

power = -1*neg_power
hydro_wrap(x_opt,true);

function negative_power = hydro_wrap(x,plotOn)

p = parameters();
p.bx = x(1);
p.bt = x(2);
p.kx = x(3);
p.kt = x(4);

power = hydro(p,plotOn);

negative_power = power;

end