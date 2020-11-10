function [x, iter, resv]= itermeth(A,b,x0,kmax,tol,P)
%ITERMETH    Un metodo iterativo generale
%  [X, ITER, RESV] = ITERMETH(A,B,X0,KMAX,TOL,P)
%  risolve iterativamente il sistema di equazioni
%  lineari A*X=B su X. La matrice A di N-per-N coef-
%  ficienti deve essere non singolare ed il termine
%  noto B deve avere lunghezza N. Se P='J' viene usato
%  il metodo di Jacobi, se P='G' viene invece selezio-
%  nato il metodo di Gauss-Seidel. Altrimenti, P e'
%  una matrice N-per-N non singolare che gioca il ruo-
%  lo di precondizionatore nel metodo del gradiente,
%  che e' un metodo di Richardson a parametro
%  dinamico. Il metodo si arresta quando il rapporto
%  fra la norma del residuo corrente e quella del
%  residuo iniziale e' minore di TOL e ITER e' il
%  numero di iterazioni effettuate. KMAX prescrive
%  il numero massimo di iterazioni consentite. Se P
%  non viene precisata, viene usato il metodo del
%  gradiente non precondizionato. In output, RESV
%  e` il vettore contenente i residui relativi
[n,n]=size(A);  resv=[ ];
if nargin == 6
 if ischar(P)==1
  if P=='J'
   L=diag(diag(A)); U=eye(n); beta=1; alpha=1;
  elseif P == 'G'
   L=tril(A); U=eye(n); beta=1; alpha=1;
  end
 else
     [L,U]=lu(P); beta = 0;
 end
else
  L = eye(n); U = L; beta = 0;
end
iter = 0;        x = x0;
r = b - A * x0;  r0 = norm(r);  res = tol+1;
while res > tol & iter < kmax
  iter = iter + 1;  z = L\r;  z = U\z;
  if beta == 0
     alpha = z'*r/(z'*A*z);
  end
  x = x + alpha*z;     r = b - A * x;
  res = norm (r) / r0; resv = [resv; res];
end
