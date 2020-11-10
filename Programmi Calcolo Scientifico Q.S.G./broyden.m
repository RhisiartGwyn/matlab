function [zero,res,niter,difv]=broyden(fun,B0,x0,...
                                 tol, kmax,varargin)
%BROYDEN Trova uno zero di un sistema di funzioni.
%   ZERO=BROYDEN(FUN,B0,X0,TOL,KMAX) approssima lo
%   zero ZERO del sistema di funzioni definite nella
%   function FUN, usando il metodo di Broyden
%   partendo da X0, dove B0 e' l'approssimazione
%   dello Jacobiano al passo 0. FUN accetta in ingresso
%   un vettore x e restituisce un vettore della stessa
%   dimensione. Se la ricerca dello zero fallisce in
%   KMAX iterazioni, il programma restituisce un mes-
%   saggio d'errore. FUN puo' essere una anonymous
%   function o una function definita in M-file.
%   ZERO=BROYDEN(FUN,B0,X0,TOL,KMAX,P1,P2,...)
%   passa i parametri P1,P2,... alla funzione
%   FUN(X,P1,P2,...).
%   [ZERO,RES,NITER,DIFV]= BROYDEN(FUN,...) restituisce
%   il valore del residuo RES in ZERO, il numero di
%   iterazioni NITER necessario per calcolare ZERO ed
%   il vettore DIFV delle norme ||x^(k+1)-x^(k)||

fx0 = fun(x0,varargin{:});
k = 0;
diff = tol+1; difv= [ ];
while diff >= tol && k < kmax
  k = k + 1;
  deltax=-B0\fx0;
  x1=x0+deltax;
  fx1=fun(x1,varargin{:});
  B0=B0+(fx1*deltax')/(deltax'*deltax);
  diff = norm(deltax);
  difv=[difv;diff];
  x0=x1; fx0=fx1;
end
zero = x1; res = norm(fx1);
if (k==kmax && diff > tol)
  fprintf(['Broyden si e'' arrestato senza aver ',...
   'soddisfatto l''accuratezza richiesta, avendo\n',...
   'raggiunto il massimo numero di iterazioni\n']);
end
niter=k;
