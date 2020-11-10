function [u]=beonestep(t,u,y,h,f,fn)
% BEONESTEP un passo del metodo di Eulero all'indietro
u = u + h*f(t,y);
