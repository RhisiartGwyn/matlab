function [alpha,niter,difv]=fixedpoint(phi,x0,tol,kmax)
%FIXEDPOINT Trova un punto fisso di una funzione.
% [ALPHA,K,DIFV]=FIXEDPOINT(PHI,X0,TOL,KMAX)
% Cerca il punto fisso ALPHA della funzione PHI
% partendo dal punto iniziale X0. TOL e KMAX sono
% rispettivamente tolleranza e numero massimo di
% iterazioni per il test d'arresto (sull'incremento).
% NITER e' il numero di iterazioni effettuate e
% DIFV e' un vettore che contiene gli errori
% |x^{k+1}-x^{k}| ad ogni passo k.
x = x0;
k = 0; err = tol+1;difv=[ ];
while err >= tol && k < kmax
   k=k+1; xnew=phi(x);
   err=abs(xnew-x); difv=[difv;err];
   x=xnew;
end
alpha = x; niter=k;
