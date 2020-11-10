function [u]=feonestep(t,y,h,f)
% FEONESTEP un passo del metodo di Eulero in avanti
u = y + h*f;
