function [zero,res,niter]=bisection(fun,a,b,tol,...
                                     kmax,varargin)
%BISECTION Trova uno zero di una funzione.
%  ZERO=BISECTION(FUN,A,B,TOL,KMAX) approssima uno
%  zero della funzione FUN nell'intervallo [A,B] con
%  il metodo di bisezione. FUN deve essere definita
%  su variabile di tipo array e puo' essere una anony-
%  mous function od una function definita in un M-file.
%  Se la ricerca dello zero di FUN fallisce, il
%  programma restituisce un messaggio d'errore.
%
%  ZERO=BISECTION(FUN,A,B,TOL,KMAX,P1,P2,...) passa
%  i parametri P1, P2,... alla funzione
%  FUN(X,P1,P2,...).
%
%  [ZERO,RES,NITER]= BISECTION(FUN,...) restituisce
%  il valore del residuo RES in ZERO ed il numero di
%  iterazioni effettuate per calcolare il valore ZERO.

x = [a, (a+b)*0.5, b];
fx = fun(x,varargin{:});
if fx(1)*fx(3) > 0
  error([' Il segno della funzione agli estremi',...
   ' dell''intervallo [A,B] deve essere diverso']);
elseif fx(1) == 0
    zero = a; res = 0; niter = 0; return
elseif fx(3) == 0
    zero = b; res = 0; niter = 0; return
end
niter = 0;
I = (b - a)*0.5;
while I >= tol && niter < kmax
   niter = niter + 1;
   if fx(1)*fx(2) <  0
      x(3) = x(2);
      x(2) = x(1)+(x(3)-x(1))*0.5;
      fx = fun(x,varargin{:});
      I = (x(3)-x(1))*0.5;
   elseif fx(2)*fx(3) < 0
      x(1) = x(2);
      x(2) = x(1)+(x(3)-x(1))*0.5;
      fx = fun(x,varargin{:});
      I = (x(3)-x(1))*0.5;
   else
       x(2) = x(find(fx==0)); I = 0;
   end
end
if  (niter==kmax && I > tol)
 fprintf(['Il metodo di bisezione si e'' arrestato',...
 ' senza soddisfare la tolleranza richiesta\n',...
 'avendo raggiunto il numero massimo di iterazioni\n']);
end
zero = x(2); x = x(2);
res = fun(x,varargin{:});
