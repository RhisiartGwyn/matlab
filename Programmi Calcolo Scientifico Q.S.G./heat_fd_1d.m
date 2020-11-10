function [xh,uh]=heat_fd_1d(xspan,tspan,nstep,mu,...
               u0,gd,f,theta,varargin)
%HEAT_FD_1D risolve l'equazione del calore con il
%  theta-metodo in tempo e differenze finite in spazio.
%  [XH,UH]=HEAT_FD_1D(XSPAN,TSPAN,NSTEP,MU,U0,GD,F,...
%  THETA) risolve l'equazione del calore
%  D U/DT - MU D^2U/DX^2 = F nel dominio
%  (XSPAN(1),XSPAN(2))x(TSPAN(1),TSPAN(2)) utilizzando
%  il theta-metodo con condizione iniziale U(X,0)=U0(X)
%  e condizioni al bordo di Dirichlet U(X,T)=GD(X,T)
%  per X=XSPAN(1) e X=XSPAN(2). MU>0 e' una costante.
%  F=F(X,T), GD=GD(X,T) e U0=U0(X) sono function handle
%  o user defined function.
%  NSTEP(1) e' il n.ro di intervalli in spazio
%  NSTEP(2) e' il n.ro di intervalli in tempo
%  XH contiene i nodi della discretizzazione
%  UH contiene la soluzione numerica al tempo TSPAN(2).
%  [XH,UH]=HEAT_FD_1D(XSPAN,TSPAN,NSTEP,MU,U0,GD,F,...
%  THETA,P1,P2,...) passa i parametri opzionali
%  P1,P2,...to alle funzioni U0,GD,F.
h  = (xspan(2)-xspan(1))/nstep(1);
dt = (tspan(2)-tspan(1))/nstep(2);
N = nstep(1)-1; e = ones(N,1);
Afd = spdiags([-e 2*e -e],[-1,0,1],N,N)/h^2;
I = speye(N);
A = I+mu*dt*theta*Afd; An = I-mu*dt*(1-theta)*Afd;
ad=theta*mu/h^2*dt; adn=(1-theta)*mu/h^2*dt;
xh = (linspace(xspan(1),xspan(2),N+2))';
xhi=xh(2:end-1);
fn = f(xhi,tspan(1),varargin{:});
un = u0(xhi,varargin{:});
R=chol(A);
gdn=gd([xspan(1);xspan(2)],tspan(1),varargin{:});
for t = tspan(1)+dt:dt:tspan(2)
    fn1 = f(xhi,t,varargin{:});
    b = An*un+dt*(theta*fn1+(1-theta)*fn);
    gdn1 = gd([xspan(1);xspan(2)],t,varargin{:});
    b([1,end]) = b([1,end])+ad*gdn1+adn*gdn;
    u = R'\b; u = R\u; fn = fn1; un = u;
    gdn=gdn1;
end
uh=[gdn1(1);u;gdn1(end)];
