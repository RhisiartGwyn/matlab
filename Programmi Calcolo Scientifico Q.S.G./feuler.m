function [t,u]=feuler(odefun,tspan,y0,Nh)
%FEULER  Risolve equazioni differenziali
%   usando il metodo di Eulero in avanti.
%   [T,Y] = FEULER(ODEFUN,TSPAN,Y0,NH) con
%   TSPAN = [T0,TF] integra il sistema di equazioni
%   differenziali y' = f(t,y) dal tempo T0 a TF con
%   condizione iniziale Y0 usando il metodo di Eulero
%   in avanti su una griglia equispaziata di NH
%   intervalli.
%   La funzione ODEFUN(T,Y) deve ritornare un vettore
%   contenente f(t,y), della stessa dimensione di y.
%   Ogni riga del vettore soluzione Y corrisponde ad
%   un istante temporale del vettore colonna T.
h=(tspan(2)-tspan(1))/Nh;
t=linspace(tspan(1),tspan(2),Nh+1)';
y0=y0(:); % genera sempre un vettore colonna
d=length(y0);
u=zeros(Nh+1,d);
u(1,:)=y0.'; % trasposta anche di variabili complesse
for n=1:Nh
    wn=u(n,:).';
    w=wn+h*odefun(t(n),wn);
    u(n+1,:)=w.';
end
