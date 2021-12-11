clear;close all;clc

tf = 1;
s0 = [0 pi/10 0 0 0]';
dt = 0.01;
p = parameters;

%fun = @dynamics;
fun = @dynamicsLagrange;
[state,time] = ode45wrap(fun,tf,s0,p,dt);

r = state(1,:);
theta = state(2,:);
power = state(5,:);
size = length(time);
%[points,CG,CB] = get_points(p,size,theta);
[points,CG,CB] = get_points(p,size,theta,r);

saveMovie = false;
animate(p,time,points,CB,CG,saveMovie);

figure;plot(time,r,time,theta,time,power)
legend('x','\theta','Power')
xlabel('time')

function dsdt = dynamicsLagrange(t, s, p)
    dsdt = [0 0 0 0 0]'; % s = [r th rdot thdot power]
    dsdt(1) = s(3);
    dsdt(2) = s(4);
    dsdt(3:4) = generated_accels(t, s', p.I,p.mu,p.m,p.kt,p.kx,p.g,...
                                        p.rho,p.L,p.W,p.h0,p.Fx,p.Ft,...
                                        p.w,p.bx,p.bt,p.nu,p.Cd );
    dsdt(5) = generated_power(t, s',    p.I,p.mu,p.m,p.kt,p.kx,p.g,...
                                        p.rho,p.L,p.W,p.h0,p.Fx,p.Ft,...
                                        p.w,p.bx,p.bt,p.nu,p.Cd );
end


function dsdt = dynamics(t, s, p)

u = sin(p.w * t);
dsdt = A_matrix(s,p) * s + B_matrix(s,p) * u;

end

function A = A_matrix(s,p)

A = [0 0 1 0;
    0 0 0 1;
    -get_kx(s,p) / get_m(s,p) 0 -get_bx(s,p) 0;
    0 -get_kt(s,p) / get_I(s,p) 0 -get_bt(s,p)];

end

function B = B_matrix(s,p)

B = [0;
    0; 
    get_Fx(s,p) / get_m(s,p); 
    get_Tt(s,p) / get_I(s,p)];

end

function kx = get_kx(s,p)
	kx = p.kx;
end

function kt = get_kt(s,p)
    kt = p.kt;
end

function m = get_m(s,p)
    m = p.m;
end

function I = get_I(s,p)
    I = p.I;
end

function bx = get_bx(s,p)
    bx = p.bx;
end

function bt = get_bt(s,p)
    bt = p.bt;
end

function Fx = get_Fx(s,p)
    Fx = p.Fx;
end

function Tt = get_Tt(s,p)
    Tt = p.Ft;
end

