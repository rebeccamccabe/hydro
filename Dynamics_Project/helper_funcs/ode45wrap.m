function [state,time] = ode45wrap(odefun,tfinal,x0,p,dt)
%wrapper around ode45 for parameter struct and evaluation of solution
%   [state,time] = ode45wrap(odefun,tfinal,x0,p,dt)

    sol = ode45(@(t,x)odefun(t,x,p), [0 tfinal], x0);
    time = 0:dt:tfinal;
    state = deval(sol,time);
    
end