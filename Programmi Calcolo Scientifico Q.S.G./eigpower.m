function [lambda,x,iter]=eigpower(A,tol,kmax,x0)
%EIGPOWER    Approssima l'autovalore di modulo massimo
%            di una matrice.
%   LAMBDA = EIGPOWER(A) calcola con il metodo delle
%   potenze l'autovalore di una matrice A di modulo
%   massimo a partire da un dato iniziale pari al
%   vettore unitario.
%   LAMBDA = EIGPOWER(A,TOL,KMAX,X0) arresta il metodo
%   quando la differenza  fra due iterate consecutive
%   e' minore di TOL (il valore di default e' 1.E-06)
%   o quando il massimo numero di iterazioni KMAX (il
%   valore di default e' 100) e' stato raggiunto.
%   [LAMBDA,X,ITER] = EIGPOWER(A,TOL,KMAX,X0)
%   restituisce anche l'autovettore unitario X tale
%   che A*X=LAMBDA*X ed il numero di iterazioni
%   effettuate per calcolare X.
[n,m] = size(A);
if n ~= m, error('Solo per matrici quadrate'); end
if nargin == 1
   tol = 1.e-06;   x0 = ones(n,1);   kmax = 100;
end
x0 = x0/norm(x0);
pro = A*x0;
lambda = x0'*pro;
err = tol*abs(lambda) + 1;
iter = 0;
while err>tol*abs(lambda) & abs(lambda)~=0 & iter<=kmax
   x = pro;              x = x/norm(x);
   pro = A*x;            lambdanew = x'*pro;
   err = abs(lambdanew - lambda);
   lambda = lambdanew;   iter = iter + 1;
end
