function [radici,iter]=newtonhorner(a,x0,tol,kmax)
%NEWTONHORNER Metodo di Newton-Horner
%   [RADICI,ITER]=NEWTONHORNER(A,X0) calcola le
%   radici del polinomio
%    P(X) = A(1)*X^N + A(2)*X^(N-1) + ...
%                                + A(N)*X + A(N+1)
%   con il metodo di Newton-Horner a partire dal dato
%   iniziale X0. Il metodo si arresta per ogni radice
%   al massimo dopo 100 iterazioni o dopo che il valore
%   assoluto della differenza fra due iterate conse-
%   cutive e' minore di 1.e-04.
%   [RADICI,ITER]=NEWTONHORNER(A,X0,TOL,KMAX) consente
%   di definire la tolleranza sul criterio d'arresto
%   ed il numero massimo di iterazioni.

if nargin == 2
    tol = 1.e-04;
    kmax = 100;
elseif nargin == 3
    kmax = 100;
end
n=length(a)-1;
radici = zeros(n,1);
iter = zeros(n,1);

for k = 1:n
    % Iterazioni di Newton
    niter = 0;
    x = x0;
    diff = tol + 1;
    while niter < kmax && diff >= tol
        [pz,b] = horner(a,x);
        [dpz,b] = horner(b,x);
        xnew = x - pz/dpz;
        diff = abs(xnew-x);
        niter = niter + 1;
        x = xnew;
    end

    if (niter==kmax && diff> tol)
      fprintf([' Il metodo non converge nel numero',...
              ' massimo di iterazioni\n ']);
    end
    % Deflazione
    [pz,a] = horner(a,x);
    radici(k) = x;
    iter(k) = niter;
end
