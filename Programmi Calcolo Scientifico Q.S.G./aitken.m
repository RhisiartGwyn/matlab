function [x,niter]=aitken(phi,x0,tol,kmax,varargin)
%AITKEN Estrapolazione di Aitken
%   [ALPHA,NITER]=AITKEN(PHI,X0) calcola un'appros-
%   simazione di un punto fisso ALPHA della funzione
%   PHI a partire dal dato iniziale X0 con il metodo
%   di estrapolazione di Aitken. Il metodo si arresta
%   dopo 100 iterazioni o dopo che il valore assoluto
%   della differenza fra due iterate consecutive e'
%   minore di 1.e-04. PHI puo' essere una anonymous
%   function, o una function definita in un M-file.
%   [ALPHA,NITER]=AITKEN(PHI,X0,TOL,KMAX) consente di
%   definire la tolleranza sul criterio d'arresto ed
%   il numero massimo di iterazioni.

if nargin == 2
    tol = 1.e-04;
    kmax = 100;
elseif nargin == 3
    kmax = 100;
end
x = x0;
diff = tol + 1;
k = 0;
while k < kmax && diff >= tol
    gx = phi(x,varargin{:});
    ggx = phi(gx,varargin{:});
    xnew = (x*ggx-gx^2)/(ggx-2*gx+x);
    diff = abs(x-xnew);
    x = xnew;
    k = k  + 1;
end
niter=k;
if (k==kmax && diff>tol)
    fprintf([' Il metodo non converge nel numero',...
            ' massimo di iterazioni\n ']);
end
