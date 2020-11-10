function [lambda,x,iter]=invshift(A,mu,tol,kmax,x0)
%INVSHIFT    Approssima l'autovalore di modulo minimo.
%   LAMBDA = INVSHIFT(A) calcola con il metodo delle
%   potenze l'autovalore di una matrice A di modulo
%   minimo a partire da un dato iniziale pari al
%   vettore unitario.
%   LAMBDA = INVSHIFT(A,MU) calcola l'autovalore di A
%   piu' vicino ad un dato numero (reale o complesso)
%   MU.
%   LAMBDA = INVSHIFT(A,MU,TOL,KMAX,X0) arresta il
%   metodo quando la differenza  fra due iterate
%   consecutive e' minore di TOL (il valore di
%   default e' 1.E-06) o quando il massimo numero di
%   iterazioni KMAX (il valore di default e' 100) e'
%   stato raggiunto.
%   [LAMBDA,X,ITER] = INVSHIFT(A,MU,TOL,KMAX,X0)
%   restituisce anche l'autovettore unitario X tale che
%   A*X=LAMBDA*X ed il numero di iterazioni effettuate
%   per calcolare X.
[n,m]=size(A);
if n ~= m, error('Solo per matrici quadrate'); end
if nargin == 1
  x0 = rand(n,1); kmax = 100; tol = 1.e-06; mu = 0;
elseif nargin == 2
  x0 = rand(n,1); kmax = 100; tol = 1.e-06;
end
[L,U]=lu(A-mu*eye(n));
if norm(x0) == 0
   x0 = rand(n,1);
end
x0=x0/norm(x0);
z0=L \ x0;
pro=U \ z0;
lambda=x0'*pro;
err=tol*abs(lambda)+1;        iter=0;
while err>tol*abs(lambda) & abs(lambda)~=0 & iter<=kmax
   x = pro; x = x/norm(x);
   z=L \ x;    pro=U \ z;
   lambdanew = x'*pro;
   err = abs(lambdanew - lambda);
   lambda = lambdanew;
   iter = iter + 1;
end
lambda = 1/lambda + mu;
