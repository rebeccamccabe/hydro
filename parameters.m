
function p = parameters()

p = struct('m', 5, 'g', 9.8, 'rho', 1000, 'Cd', .01, ...
    'L', 1, 'H', 5, 'h0', 2, 'W', 1);
p.I = 1/12 * p.m * (p.L^2 + p.H^2);

p.mu = 1;
p.kx = 1;
p.kt = 1;
p.Fx = 10;
p.Ft = 10;
p.bx = 1;
p.bt = 1;
p.nu = 1;

p.Hs = .5;
p.w = 5;

assert(p.h0 > p.Hs,'Error: bottom of wave is below bottom of device')
assert(p.H - p.h0 > p.Hs,'Error: top of wave is above top of device')
assert(p.rho * p.W * p.L * p.h0 > p.m,'Error: WEC is not positively buoyant')

end