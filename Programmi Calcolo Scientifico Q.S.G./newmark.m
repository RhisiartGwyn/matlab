function [t,u]=newmark(odefun,tspan,y0,Nh,param,...
                       varargin)
%NEWMARK  Risolve equazioni differenziali del II ord.
%  con il metodo di Newmark.
%  [T,Y]=NEWMARK(ODEFUN,TSPAN,Y0,NH,PARAM) con
%  TSPAN=[T0 TF] integra il sistema di equazioni dif-
%  ferenziali y''= f(t,y,y') dal tempo T0 a TF con
%  condizioni iniziali Y0=(y(t0),y'(t0)) utilizzando
%  il metodo di Newmark su una griglia equispaziata di
%  NH intervalli. Il vettore PARAM contiene, in ordine,
%  i parametri zeta e theta del metodo di Newmark.
%  La funzione ODEFUN(T,Y) deve ritornare uno scalare,
%  la variabile Y e' un array che contiene la funzione
%  soluzione in prima componente e la sua derivata in
%  seconda componente.
%  Ogni riga del vettore soluzione Y corrisponde ad
%  un istante temporale del vettore colonna T.
%  Richiama broyden.m per la risoluzione del sistema
%  non lineare nel caso il metodo sia implicito.
%  Se sono assegnati solo 4 parametri
%  di input, tolleranza e numero max di iterazioni per
%  il metodo di Broyden sono posti pari a tol=1.e-8
%  e kmax=20 rispettivamente, altrimenti
%  [T,Y] = NEWMARK(ODEFUN,TSPAN,Y0,NH,PARAM,TOL,KMAX)
%  permette di specificare anche i parametri per
%  la function broyden.m
if nargin == 5
    tol=1.e-8;  kmax=20;
else
    tol=varargin{1}; kmax=varargin{2};
end
h=(tspan(2)-tspan(1))/Nh; h2=h^2;
t=linspace(tspan(1),tspan(2),Nh+1)';
y0=y0(:); % genera sempre un vettore colonna
u=zeros(Nh+1,2);
u(1,:)=y0.'; % trasposta anche di variabili complesse
B0=eye(2); % matrice iniziale per innestare Broyden
zeta=param(1);  zeta12=0.5-zeta;
theta=param(2); theta1=1-theta;
for n=1:Nh
  wn=u(n,:).';  fn=odefun(t(n),wn);
  if theta==0 && zeta==0
    w=[wn(1)+h*wn(2)+h2*zeta12*fn; wn(2)+h*fn];
  else
    F=@(x)[x(1)-wn(1)-h*wn(2)-h2*...
           (zeta*odefun(t(n+1),x)+zeta12*fn);
           x(2)-wn(2)-h*(theta*odefun(t(n+1),x)+...
            theta1*fn)];
    w=broyden(F,B0,wn,tol,kmax);
  end
  u(n+1,:)=w.';
end
