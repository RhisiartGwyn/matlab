function [u]=cnonestep(t,u,y,h,f,fn)
% CNONESTEP un passo del metodo di Crank-Nicolson
u = u + 0.5*h*(f(t,y)+fn);
