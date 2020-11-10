function[JSf,nodes]=simpadpt(fun,a,b,tol,hmin,varargin)
%SIMPADPT   Formula adattiva di Simpson.
%   JSF = SIMPADPT(FUN,A,B,TOL,HMIN) approssima
%   l'integrale di FUN nell'intervallo (A,B)
%   garantendo che il valore assoluto dell'errore sia
%   inferiore a TOL ed utilizzando un passo variabile
%   H >= HMIN. FUN e' una function che riceve
%   in ingresso un vettore x e restituisce un vettore
%   reale.  FUN puo' essere una anonymous function o
%   una function definita in un M-file.
%   JSF = SIMPADPT(FUN,A,B,TOL,HMIN,P1,P2,...) passa
%   alla function FUN i parametri opzionali P1,P2,...
%   come FUN(X,P1,P2,...).
%   [JSF,NODES] = SIMPADPT(...) restituisce la distri-
%   buzione di nodi usati nel processo di quadratura.
A=[a,b]; N=[]; S=[]; JSf = 0; ba = 2*(b - a); nodes=[];
while ~isempty(A),
  [deltaI,ISc]=caldeltai(A,fun,varargin{:});
  if abs(deltaI) < 15*tol*(A(2)-A(1))/ba;
     JSf = JSf + ISc;    S = union(S,A);
     nodes = [nodes, A(1) (A(1)+A(2))*0.5 A(2)];
     S = [S(1), S(end)]; A = N; N = [];
  elseif A(2)-A(1) < hmin
     JSf=JSf+ISc;        S = union(S,A);
     S = [S(1), S(end)]; A=N; N=[];
     warning('Passo di integrazione troppo piccolo');
  else
     Am = (A(1)+A(2))*0.5;
     A = [A(1) Am];
     N = [Am, b];
  end
end
nodes=unique(nodes);
return

function [deltaI,ISc]=caldeltai(A,fun,varargin)
L=A(2)-A(1);
t=[0; 0.25; 0.5; 0.75; 1];
x=L*t+A(1);
L=L/6;
w=[1; 4; 1]; wp=[1;4;2;4;1];
fx=fun(x,varargin{:}).*ones(5,1);
IS=L*sum(fx([1 3 5]).*w);
ISc=0.5*L*sum(fx.*wp);
deltaI=IS-ISc;
return
