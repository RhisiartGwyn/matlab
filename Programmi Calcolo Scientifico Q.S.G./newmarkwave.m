function [xh,uh]=newmarkwave(xspan,tspan,nstep,param,...
               c,u0,v0,gd,f,varargin)
%NEWMARKWAVE risolve l'equazione delle onde
%       con il metodo di Newmark
%   [XH,UH]=NEWMARKWAVE(XSPAN,TSPAN,NSTEP,PARAM,C,...
%   U0,V0,GD,F)
%   risolve l'equazione delle onde
%      D^2 U/DT^2 - C D^2U/DX^2 = F in
%   (XSPAN(1),XSPAN(2)) X (TSPAN(1),TSPAN(2)) con il
%   metodo di Newmark, con condizioni iniziali
%   U(X,0)=U0(X), DU/DX(X,0)=V0(X) e condizioni al
%   bordo di Dirichlet U(X,T)=GD(X,T) in X=XSPAN(1)
%   ed in X=XSPAN(2). C e' una costante positiva.
%   NSTEP(1) e' il numero di intervalli di spazio.
%   NSTEP(2) e' il numero di intervalli in tempo.
%   PARAM(1)=ZETA e PARAM(2)=THETA.
%   U0(X), V0(X), GD(X,T) e F(X,T) possono essere defi-
%   nite come function handle o user defined function.
%   XH contiene i nodi della discretizzazione in spazio
%   UH contiene la soluzione numerica al tempo TSPAN(2).
%   [XH,UH]=NEWMARKWAVE(XSPAN,TSPAN,NSTEP,PARAM,C,...
%   U0,V0,GD,F,P1,P2,...) passa i parametri addizio-
%   nali P1,P2,... alle funzioni U0,V0,GD,F.
h  = (xspan(2)-xspan(1))/nstep(1);
dt = (tspan(2)-tspan(1))/nstep(2);
zeta = param(1);  theta = param(2);
N = nstep(1)+1; e = ones(N,1);
Afd = spdiags([e -2*e e],[-1,0,1],N,N)/h^2;
I = speye(N);
A = I-c*dt^2*zeta*Afd;
An = I+c*dt^2*(0.5-zeta)*Afd;
A(1,:) = 0; A(1,1) = 1;
A(N,:) = 0; A(N,N) = 1;
xh = (linspace(xspan(1),xspan(2),N))';
fn = f(xh,tspan(1),varargin{:});
un = u0(xh,varargin{:});
vn = v0(xh,varargin{:});
[L,U]=lu(A);
alpha = dt^2*zeta;
beta = dt^2*(0.5-zeta);
theta1 = 1-theta;
for t = tspan(1)+dt:dt:tspan(2)
    fn1 = f(xh,t,varargin{:});
    rhs = An*un+dt*I*vn+alpha*fn1+beta*fn;
    temp = gd([xspan(1),xspan(2)],t,varargin{:});
    rhs([1,N]) = temp;
    uh = L\rhs;    uh = U\uh;
    v = vn + dt*((1-theta)*(c*Afd*un+fn)+...
        theta*(c*Afd*uh+fn1));
    fn = fn1;    un = uh;    vn = v;
end
