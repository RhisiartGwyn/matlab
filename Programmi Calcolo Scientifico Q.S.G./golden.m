function [xmin,fmin,iter]=golden(fun,a,b,tol,...
                                     kmax,varargin)
%GOLDEN Trova un punto di minimo di una funzione.
%  XMIN=GOLDEN(FUN,A,B,TOL,KMAX) approssima un
%  punto di minimo della funzione FUN nell'intervallo
%  [A,B] con il metodo della sezione aurea.
%  Se la ricerca del punto di minimo di FUN fallisce,
%  il programma restituisce un messaggio d'errore.
%  FUN puo' essere una anonymous function
%  od una function definita in un M-file.
%  XMIN=GOLDEN(FUN,A,B,TOL,KMAX,P1,P2,...) passa
%  i parametri P1, P2,... alla funzione
%  FUN(X,P1,P2,...).
%  [XMIN,FMIN,ITER]= GOLDEN(FUN,...) restituisce
%  il valore della funzione FUN in XMIN e il numero
%  di iterazioni effettuate per calcolare XMIN.
phi=(1+sqrt(5))/2; phi1=1/phi; phi2=1/(phi+1);
c=a+phi2*(b-a); d=a+phi1*(b-a);
err=tol+1; k=0;
while err>tol && k< kmax
  if(fun(c) >= fun(d))
    a=c; c=d; d=a+phi1*(b-a);
  else
    b=d; d=c; c=a+phi2*(b-a);
  end
  k=k+1; err=abs(b-a)/(abs(c)+abs(d));
end
xmin=(a+b)/2; fmin=fun(xmin); iter=k;
if  (iter==kmax && err > tol)
 fprintf(['Il metodo della sez. aurea si e'' \n',...
 'arrestato senza soddisfare la tolleranza \n',...
 'richiesta avendo raggiunto il numero massimo \n',...
 ' di iterazioni\n']);
end
