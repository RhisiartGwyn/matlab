function [zero,res,niter,difv]=newton(fun,dfun,x0,tol,...
                                 kmax,varargin)
%NEWTON Trova uno zero di una funzione.
%   ZERO=NEWTON(FUN,DFUN,X0,TOL,KMAX) approssima lo
%   zero ZERO della funzione definita nella function
%   FUN, continua e derivabile, usando il metodo di
%   Newton e partendo da X0. Se la ricerca
%   dello zero fallisce in KMAX iterazioni, il pro-
%   gramma restituisce un messaggio d'errore.
%   FUN e DFUN possono essere anonymous
%   function o function definite in M-file.
%   ZERO=NEWTON(FUN,DFUN,X0,TOL,KMAX,P1,P2,...) passa
%   i parametri P1,P2,... alle funzioni
%   FUN(X,P1,P2,...) e DFUN(X,P1,P2,...).
%   [ZERO,RES,NITER,DIFV]= NEWTON(FUN,...) restituisce
%   il valore del residuo RES in ZERO, il numero di
%   iterazioni NITER necessario per calcolare ZERO  ed
%   il vettore DIFV degli incrementi |x^(k+1)-x^k|

x = x0;
fx = fun(x,varargin{:});
dfx = dfun(x,varargin{:});
k = 0; diff = tol+1; difv=[ ];
while diff >= tol && k < kmax
   k = k + 1;
   diff = - fx/dfx;
   x = x + diff;
   diff = abs(diff); difv=[difv; diff];
   fx = fun(x,varargin{:});
   dfx = dfun(x,varargin{:});
end
if (k==kmax && diff > tol)
  fprintf(['Newton si e'' arrestato senza aver ',...
   'soddisfatto l''accuratezza richiesta, avendo\n',...
   'raggiunto il massimo numero di iterazioni\n']);
end
zero = x; res = fx; niter=k;
