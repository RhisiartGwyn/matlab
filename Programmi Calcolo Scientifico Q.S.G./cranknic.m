function [t,u]=cranknic(odefun,tspan,y0,Nh,varargin)
%CRANKNIC  Risolve equazioni differenziali
%  usando il metodo di Crank-Nicolson.
%  [T,Y]=CRANKNIC(ODEFUN,TSPAN,Y0,NH) con
%   TSPAN = [T0,TF] integra il sistema di equazioni
%   differenziali y' = f(t,y) dal tempo T0 a TF con
%   condizione iniziale Y0 usando il metodo di
%   Crank-Nicolson su una griglia equispaziata di NH
%   intervalli.
%   La funzione ODEFUN(T,Y) deve ritornare un vettore
%   contenente f(t,y), della stessa dimensione di y.
%   Ogni riga del vettore soluzione Y corrisponde ad
%   un istante temporale del vettore colonna T.
%   Richiama broyden.m per la risoluzione del sistema
%   non lineare. Se sono assegnati solo 4 parametri
%   di input, tolleranza e numero max di iterazioni per
%   il metodo di Broyden sono posti pari a tol=1.e-8
%   e kmax=20 rispettivamente, altrimenti
%   [T,Y] = CRANKNIC(ODEFUN,TSPAN,Y0,NH,TOL,KMAX)
%   permette di specificare anche i parametri per
%   function broyden.m.
if nargin == 4
    tol=1.e-8;  kmax=20;
else
    tol=varargin{1}; kmax=varargin{2};
end
h=(tspan(2)-tspan(1))/Nh;  h2=h/2;
t=linspace(tspan(1),tspan(2),Nh+1)';
y0=y0(:); % genera sempre un vettore colonna
d=length(y0);
u=zeros(Nh+1,d);
B0=eye(d); % matrice iniziale per innestare Broyden
u(1,:)=y0.'; % trasposta anche di variabili complesse
for n=1:Nh
    wn=u(n,:).';
    F=@(x)x-wn-h2*(odefun(t(n),wn)+odefun(t(n+1),x));
    w=broyden(F,B0,wn,tol,kmax);
    u(n+1,:)=w.';
end
