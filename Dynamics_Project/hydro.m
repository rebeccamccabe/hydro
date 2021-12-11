function P = hydro(p,plotOn)

if nargin==0
    close all;clc
    p = parameters();
    plotOn = true;
end

tf = 5;
s0 = [0 0 0 0 0]';
dt = 0.01;

%fun = @dynamicsNewtonian;
fun = @dynamicsLagrange;
[state,time] = ode45wrap(fun,tf,s0,p,dt);

power = state(5,:);
P = mean(power);

if plotOn
    r = state(1,:);
    theta = state(2,:);
    size = length(time);
    [points,CG,CB] = get_points(p,size,theta,r);

    saveMovie = false;
    plot_title = ['WEC response to ', num2str(p.Hs), ' m wave'];
    animate(p,time,points,CB,CG,saveMovie,plot_title);

    figure;plot(time,r,time,theta,time,power)
    legend('r','\theta','Power')
    xlabel('time')
end

end

function dsdt = dynamicsLagrange(t, s, p)
    dsdt = [0 0 0 0 0]'; % s = [r th rdot thdot power]
    dsdt(1) = s(3);
    dsdt(2) = s(4);
    dsdt(3:4) = generated_accels(t, s', p.mu,p.m,p.kt,p.kx,p.g,...
                                        p.rho,p.L,p.W,p.h0,p.Fx,p.Ft,...
                                        p.w,p.bx,p.bt,p.nu,p.Cd );
    dsdt(5) = generated_power(t, s',    p.mu,p.m,p.kt,p.kx,p.g,...
                                        p.rho,p.L,p.W,p.h0,p.Fx,p.Ft,...
                                        p.w,p.bx,p.bt,p.nu,p.Cd );
end


function dsdt = dynamicsNewtonian(t, s, p)

u = sin(p.w * t);
A = A_matrix(s,p);
B = B_matrix(s,p);

% account for s=0 causing infinite stiffness/inertia
bad_idx = abs(A) == Inf;
A(bad_idx) = 1e6 * sign(A(bad_idx)); 
bad_idx = abs(B) == Inf;
B(bad_idx) = 1e6 * sign(B(bad_idx)); 

dsdt = A * s + B * u;
assert(~any(isnan(dsdt)),'NaN error')

end

function A = A_matrix(s,p)

A = [0 0 1 0 0;
    0 0 0 1 0;
    -get_Kx(s,p) / get_m(s,p) 0 -get_Bx(s,p) 0 0;
    0 -get_Kt(s,p) / get_I(s,p) 0 -get_Bt(s,p) 0];

A(5,:) = p.kx*A(1,:) + p.kt*A(2,:) + p.bx*A(3,:) + p.bt*A(4,:);

end

function B = B_matrix(s,p)

B = [0;
    0; 
    get_Fx(s,p) / get_m(s,p); 
    get_Tt(s,p) / get_I(s,p); 
    0];

end

function kx = get_Kx(s,p)
    r = s(1) + p.H/2;
    theta = s(2);
    if theta == 0
        kx = p.kx;
    else  
        kx = p.kx + p.rho*p.g*p.W*p.L*p.h0* (1-cos(theta))/(r*cos(theta));
    end
    assert(kx > 0,'Error: negative stiffness')
end

function kt = get_Kt(s,p)
    r = s(1) + p.H/2;
    theta = s(2);
    if theta == 0
        kt = p.kt; % account for divide by zero case
    else
        kt = p.kt - sin(theta)/theta * p.g/2 * (2*p.m*r - p.rho*p.L*p.W*p.h0^2/cos(theta)^2);
    end
    assert(kt > 0,'Error: negative stiffness')
end

function m = get_m(s,p)
    m = p.m;
end

function I = get_I(s,p)
    r = s(1) + p.H/2; % +H/2 is required because my animation plots assuming r equilibrium = 0 instead of H/2
    I = 4/3 * p.m * (r)^2; 
end

function bx = get_Bx(s,p)
    r_dot = s(3);
    bx = p.bx + 1/2 * p.rho * p.Cd * p.W * p.L * abs(r_dot);
end

function bt = get_Bt(s,p)
    theta_dot = s(4);
    bt = (p.bt + p.nu) + 1/2 * p.rho * p.Cd * p.W * p.h0^4 * abs(theta_dot);
end

function Fx = get_Fx(s,p)
    Fx = p.Fx;
end

function Tt = get_Tt(s,p)
    Tt = p.Ft;
end

