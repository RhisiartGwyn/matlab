function [x,relres,k,res]=gradient(A,b,tol,kmax,RP,x0)
% GRADIENT metodo del gradiente (precondizionato).
%  [X,RELRES,K,RES]=GRADIENT(A,B,TOL,KMAX,RP,X0)
%  risolve il sistema A*X=B con il metodo del gradiente
%  (precondizionato). X0 e` il vettore iniziale,
%  TOL specifica la tolleranza per il test d'arresto.
%  KMAX indica il numero massimo di iterazioni ammesse.
%  RP contiene una matrice triangolare superiore t.c.
%  il precondizionatore e` definito come P = RP'*RP.
%  [X,RELRES,K,RES]=GRADIENT(A,B,TOL,KMAX,[ ],X0) ri-
%  chiama il metodo del gradiente non precondizionato.
%  RELRES e` la norma del residuo corrispondende alla
%  soluzione calcolata diviso la norma del vettore B.
%  K e` il numero di iterazioni effettuate dal metodo,
%  RES e` un vettore contenente la norma dei residui
%  ad ogni passo del metodo

k = 0;
bnorm = norm(b);
res=[ ];
if bnorm==0
x=zeros(n,1); relres=0;
else
r=b-A*x0; x=x0;
% sistema sul precondizionatore
if ~isempty(RP); z=RP' \ (RP\r); else  z=r; end
zr=z'*r; res0=sqrt(zr); relres=zr;
while relres > tol && k< kmax
    v=A*z;
    alpha=zr/(z'*v);
    x=x+alpha*z;
    r=r-alpha*v;
    if ~isempty(RP); z=RP' \ (RP\r); else  z=r; end
    zr=z'*r;
    relres=sqrt(zr)/res0;
    res=[res; relres];
    k=k+1;
end
end
