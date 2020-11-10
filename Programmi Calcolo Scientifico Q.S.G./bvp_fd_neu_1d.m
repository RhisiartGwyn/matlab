function [xh,uh]=bvp_fd_neu_1d(a,b,N,mu,eta,sigma,...
                    bvpfun,gamma,delta,varargin)
% BVP_FD_NEU_1D Risolve il problema di diffusione-tra-
% sporto-reazione e condizioni di Neumann su (a,b)
% con differenze finite centrate e N+2 nodi equispa-
% ziati.
% [xh,uh]=bvp_fd_neu_1d(a,b,N,mu,eta,sigma,...
%                   bvpfun,gamma,delta)
%
h =(b-a)/(N+1); xh = (linspace(a,b,N+2))';
hm =mu/h^2;  hd = eta/(2*h);
e =ones(N+2,1);
Afd =spdiags([-(hm+hd)*e (2*hm+sigma)*e (-hm+hd)*e],...
               -1:1, N, N);
Afd(1,1)=3/(2*h); Afd(1,2)=-2/h; Afd(1,3)=1/(2*h);
Afd(N+2,N+2)=3/(2*h); Afd(N+2,N+1)=-2/h;
Afd(N+2,N)=1/(2*h);
f =bvpfun(xh,varargin{:});
f(1)=gamma; f(N+2)=delta;
uh = Afd\f;
