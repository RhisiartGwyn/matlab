function [xh,uh]=bvp_fe_dir_1d(a,b,N,mu,eta,sigma,...
                    bvpfun,ua,ub,varargin)
%BVP_FE_DIR_1D Risolve un problema ai limiti
%   [XH,UH]=BVP_FE_DIR_1D(A,B,N,MU,ETA,SIGMA,BVPFUN,...
%           UA,UB) con il metodo degli elementi
%   finiti di grado 1 con passo h=1/N il problema
%   -MU*D(DU/DX)/DX + ETA*DU/DX + SIGMA*U = BVPFUN
%   sull'intervallo (A,B) con condizioni al bordo
%   U(A)=UA e U(B)=UB. Il termine di trasporto e'
%   trattato con schema upwind.
%   BVPFUN puo' essere un function handle o una
%   user defined function.
%   In output, XH e UH contengono, risp., i nodi e la
%   soluzione numerica, inclusi i valori al bordo.
h = (b-a)/(N+1); xh = (linspace(a,b,N+2))';
Pe=eta*h/(2*mu); hm = mu*(1+Pe)/h; hd = eta/2;
hs=sigma*h/6; e =ones(N,1);
% stiffness matrix
Afe =spdiags([(-hm-hd+hs)*e (2*hm+4*hs)*e...
             (-hm+hd+hs)*e], -1:1, N, N);
f=quad_mp(bvpfun,xh);
f(1) =  f(1)+mu*ua/h;
f(end) = f(end)+mu*ub/h;
uh = Afe\f;
uh=[ua; uh; ub];
return
function [f]=quad_mp(bvpfun,xh);
N=length(xh)-2; h=xh(2)-xh(1);
f=zeros(N,1);
for j=1:N
    f(j)=h/2*(bvpfun((xh(j)+xh(j+1))/2)+...
        bvpfun((xh(j+1)+xh(j+2))/2));
end
return
