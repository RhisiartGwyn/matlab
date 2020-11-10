function [tt,u]=rk2(odefun,tspan,y0,Nh)
% RK2 Risolve un sistema di e.d.o. con Runge-Kutta2
% (noto anche come metodo di Heun) e passo costante
h=(tspan(2)-tspan(1))/Nh;  hh=h*0.5;
tt=linspace(tspan(1),tspan(2),Nh+1)';
y0=y0(:); % genera sempre un vettore colonna
d=length(y0);
u=zeros(Nh+1,d);
u(1,:)=y0.'; % trasposta anche di variabili complesse
for n=1:Nh
    wn=u(n,:).';
    k1=odefun(tt(n),wn);
    y=wn+h*k1;
    k2=odefun(tt(n+1),y);
    w=wn+hh*(k1+k2);
    u(n+1,:)=w.';
end
