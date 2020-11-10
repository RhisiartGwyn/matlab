function [x,err,iter]= gaussnewton(r,jr,x0,tol,...
    kmax,varargin)
%GAUSSNEWTON  Minimi quadrati non lineari
%  [X,ERR,ITER]=GAUSSNEWTON(R,JR,X0,TOL,KMAX)
%  Risolve il problema dei minimi quadrati non lineari
%  mediante il metodo di Gauss-Newton.
%  R e JR sono function handle associati alla funzione
%  R ed allo Jacobiano di R, in cui la prima variabile
%  in input e` X, di seguito possono essere dati
%  parametri opzionali. X0 e' il punto iniziale
%  della successione. TOL e' la tolleranza per il
%  test d'arresto e KMAX e' il numero massimo di
%  iterazioni.
err=tol+1; k=0; xk=x0(:);
rk=r(xk,varargin{:}); jrk=jr(xk,varargin{:});
while err>tol && k< kmax
[Q,R]=qr(jrk,0); dk=-R \ (Q'*rk);
xk1=xk+dk;
rk1=r(xk1,varargin{:});
jrk1=jr(xk1,varargin{:});
k=k+1;  err=norm(xk1-xk);
xk=xk1; rk=rk1; jrk=jrk1;
end
x=xk; iter=k;
if (k==kmax && err > tol)
fprintf(['GaussNewton si e'' arrestato senza aver ',...
   'soddisfatto l''accuratezza richiesta, avendo\n',...
   'raggiunto il massimo numero di iterazioni\n']);
end
