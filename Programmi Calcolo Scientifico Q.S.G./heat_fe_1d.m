function [xh,uh]=heat_fe_1d(xspan,tspan,nstep,mu,...
               u0,gd,f,theta,varargin)
%HEAT_FE_1D risolve l'equazione del calore con il
%  theta-metodo in tempo ed elementi finiti in spazio.
%  [XH,UH]=HEAT_FE_1D(XSPAN,TSPAN,NSTEP,MU,U0,GD,F,...
%          THETA) risolve l'equazione del calore
%  D U/DT - MU D^2U/DX^2 = F nel dominio
%  (XSPAN(1),XSPAN(2))x(TSPAN(1),TSPAN(2)) utilizzando
%  il theta-metodo con condizione iniziale U(X,0)=U0(X)
%  e condizioni al bordo di Dirichlet U(X,T)=GD(X,T)
%  per X=XSPAN(1) e X=XSPAN(2). MU e' una costante
%  positiva. F=F(X,T), GD=GD(X,T) e U0=U0(X) sono
%  function handle o user defined function.
%  NSTEP(1) e' il n.ro di intervalli in spazio
%  NSTEP(2) e' il n.ro di intervalli in tempo
%  XH contiene i nodi della discretizzazione
%  UH contiene la soluzione numerica al tempo TSPAN(2).
%  [XH,UH]=HEAT_FE_1D(XSPAN,TSPAN,NSTEP,MU,U0,GD,F,...
%  THETA,P1,P2,...) passa i parametri opzionali
%  P1,P2,...to alle funzioni U0,GD,F.
h  = (xspan(2)-xspan(1))/nstep(1);
dt = (tspan(2)-tspan(1))/nstep(2);
N = nstep(1)-1; theta1=1-theta;
e = ones(N,1); hm=mu/h;
Afe =spdiags([-hm*e 2*hm*e -hm*e],-1:1,N,N);
M= spdiags([e 4*e e],-1:1, N, N)/6*h;
A = M/dt+theta*Afe; An = M/dt-theta1*Afe;
ad=-h/(6*dt)+theta*hm;  adn=h/(6*dt)+theta1*hm;
xh = (linspace(xspan(1),xspan(2),N+2))';
xhi=xh(2:end-1);
un = u0(xhi,varargin{:});
fn=quad_mpt(f,xh,tspan(1));
R=chol(A);
gdn=gd([xspan(1);xspan(2)],tspan(1),varargin{:});
for t = tspan(1)+dt:dt:tspan(2)
    fn1=quad_mpt(f,xh,t);
    b = An*un+theta*fn1+theta1*fn;
    gdn1 = gd([xspan(1);xspan(2)],t,varargin{:});
    b([1,end]) = b([1,end])+ad*gdn1+adn*gdn;
    u = R'\b; u = R\u; fn = fn1; un = u;
    gdn=gdn1;
end
uh=[gdn1(1);u;gdn1(end)];
return
function [fn]=quad_mpt(f,xh,t)
N=length(xh)-2; h=xh(2)-xh(1);
fn=zeros(N,1);
for j=1:N
    fn(j)=h/2*(f((xh(j)+xh(j+1))/2,t)+...
        f((xh(j+1)+xh(j+2))/2,t));
end
return
