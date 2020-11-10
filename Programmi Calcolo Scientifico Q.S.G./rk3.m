function [tt,u]=rk3(odefun,tspan,y0,Nh)
% RK3 Risolve un sistema di e.d.o. con Runge-Kutta3
% con passo costante
h=(tspan(2)-tspan(1))/Nh;  hh=h*0.5; h6=h/6; h2=h*2;
tt=linspace(tspan(1),tspan(2),Nh+1)';
y0=y0(:); % genera sempre un vettore colonna
d=length(y0);
u=zeros(Nh+1,d);
u(1,:)=y0.'; % trasposta anche di variabili complesse
for n=1:Nh
    wn=u(n,:).';
    k1=odefun(tt(n),wn);
    t1=tt(n)+hh; y=wn+hh*k1;
    k2=odefun(t1,y);
    y=wn+h*(2*k2-k1);
    k3=odefun(tt(n+1),y);
    w=wn+h6*(k1+4*k2+k3);
    u(n+1,:)=w.';
end
