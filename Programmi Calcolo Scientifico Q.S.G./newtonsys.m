function [x,res,niter,difv] = newtonsys(Ffun,Jfun,...
                                x0,tol, kmax, varargin)
%NEWTONSYS calcola una radice di un sistema non lineare
%  [ZERO,RES,NITER,DIFV]=NEWTONSYS(FFUN,JFUN,X0,...
%                                  TOL, KMAX)
%  calcola il vettore ZERO, radice di un sistema non
%  lineare definito nella function FFUN con matrice
%  Jacobiana definita nella function JFUN a partire
%  dal vettore X0. RES contiene il valore del residuo
%  in ZERO e NITER il numero di iterazioni necessarie
%  per calcolare ZERO. Il vettore DIFV contiene le
%  norme ||x^(k+1)-x^(k)||. FFUN e JFUN possono essere
%  anonymous functions o user-defined functions.

k = 0;
err = tol + 1; difv=[ ];
x = x0;
while err >= tol && k < kmax
    J = Jfun(x,varargin{:});
    F = Ffun(x,varargin{:});
    delta = - J\F;
    x = x + delta;
    err = norm(delta); difv=[difv; err];
    k = k + 1;
end
res = norm(Ffun(x,varargin{:}));
if (k==kmax && err> tol)
    fprintf(['Il metodo non converge nel massimo ',...
       'numero di iterazioni. L''ultima iterata\n',...
       'calcolata ha residuo relativo pari a %e\n'],F);
end
niter=k;
