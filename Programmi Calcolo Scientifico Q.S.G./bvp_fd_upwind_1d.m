function [xh,uh]=bvp_fd_upwind_1d(a,b,N,mu,eta,sigma,...
                    bvpfun,ua,ub,varargin)
% BVP_FD_UPWIND_1D Risolve il problema di diffusione-
% trasporto e condizioni di Dirichlet su (a,b)
% implementando lo schema upwind (differenze finite
% all'indietro per il termine del trasporto, con eta>0).
% [xh,uh]=bvp_fd_upwind_1d(a,b,N,mu,eta,sigma,...
%                   bvpfun,ua,ub)
%
h = (b-a)/(N+1);
xh = (linspace(a,b,N+2))';
hm = mu/h^2;
hd = eta/h;
e =ones(N,1);
Afd = spdiags([-(hm+hd)*e (2*hm+hd+sigma)*e -hm*e],...
    -1:1, N, N);
xi = xh(2:end-1);
f =bvpfun(xi,varargin{:});
f(1) =  f(1)+ua*(hm+hd);
f(end) = f(end)+ub*(hm-hd);
uh = Afd\f;
uh=[ua; uh; ub];
