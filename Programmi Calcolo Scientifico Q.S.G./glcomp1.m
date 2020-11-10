function [intgl]=glcomp1(fun,a,b,M,varargin)
% GLCOMP1 calcola un integrale con formule Gaussiane
% INTGL=GLCOMP1(FUN,A,B,M,VARARGIN)
% INTGL e' l'approssimazione dell'integrale di FUN
% sull'intervallo (A,B) calcolata con la formula com-
% posita di Gauss-Legendre di grado 1 su M intervallini
y = [-1/sqrt(3),1/sqrt(3)];
H2 = (b-a)/(2*M);
z = [a:2*H2:b];
zM = (z(1:end-1)+z(2:end))*0.5;
x = [zM+H2*y(1), zM+H2*y(2)];
f = fun(x,varargin{:});
intgl = H2*sum(f);
return
