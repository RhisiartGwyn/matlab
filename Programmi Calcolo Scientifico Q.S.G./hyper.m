function [xh,th,uh]=hyper(xspan,tspan,u0,ul,...
                          scheme,cfl,h,deltat)
% HYPER risolve un'eqz scalare iperbolica, a>0
% [XH,TH,UH]=HYPER(XSPAN,TSPAN,U0,UL,SCHEME,CFL,...
%                 H,DELTAT)
% risolve l'equazione differenziale iperbolica scalare
%       DU/DT+ A * DU/DX=0
% in (XSPAN(1),XSPAN(2))x(TSPAN(1),TSPAN(2))
% con condizione iniziale  U(X,0)=U0(X) e
% condizione al bordo U(T)=UL(T) assegnata in XSPAN(1)
% con vari schemi alle differenze finite.
% scheme = 1 Lax - Friedrichs
%          2 Lax - Wendroff
%          3 Upwind
% La velocita' di propagazione A non e' richiesta
% esplicitamente, essendo CFL = A * DELTAT / DELTAX
% In output XH e' il vettore della discretizzazione
% in x; TH e' il vettore della discretizzazione in t
% UH e' una matrice che contiene la soluzione numerica:
% UH(n,:) contiene la sol all'istante temporale TT(n)
% U0 e UL possono essere inline o anonymous function o
% function definite tramite M-file.

Nt=(tspan(2)-tspan(1))/deltat+1;
th=linspace(tspan(1),tspan(2),Nt);
Nx=(xspan(2)-xspan(1))/h+1;
xh=linspace(xspan(1),xspan(2),Nx);
u=zeros(Nt,Nx); cfl2=cfl*0.5; cfl21=1-cfl^2;
cflp1=cfl+1; cflm1=cfl-1;
uh(1,:)=u0(xh);
for n=1:Nt-1
 uh(n+1,1)=ul(th(n+1));
 if scheme == 1
% Lax Friedrichs
    for j=2:Nx-1
      uh(n+1,j)=0.5*(-cflm1*uh(n,j+1)+cflp1*uh(n,j-1));
    end
    j=Nx;
    uh(n+1,j)=0.5*(-cflm1*(2*uh(n,j)-uh(n,j-1))+...
        cflp1*uh(n,j-1));
 elseif scheme == 2
% Lax Wendroff
    for j=2:Nx-1
     uh(n+1,j)=cfl21*uh(n,j)+...
         cfl2*(cflm1*uh(n,j+1)+cflp1*uh(n,j-1));
    end
    j=Nx;
    uh(n+1,j)=cfl21*uh(n,j)+...
     cfl2*(cflm1*(2*uh(n,j)-uh(n,j-1))+cflp1*uh(n,j-1));
 elseif scheme ==3
% Upwind
   for j=2:Nx
        uh(n+1,j)=-cflm1*uh(n,j)+cfl*uh(n,j-1);
    end
 end
end
