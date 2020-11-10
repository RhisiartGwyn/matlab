function [xh,uh]=bvp_fd_dir_1d(a,b,N,mu,eta,sigma,...
                    bvpfun,ua,ub,varargin)
%BVP_FD_DIR_1D Risolve un problema ai limiti
%  [XH,UH]=BVP_FD_DIR_1D(A,B,N,MU,ETA,SIGMA,BVPFUN,...
%          UA,UB) risolve con il metodo delle
%  differenze finite centrate in N nodi equispaziati
%  interni ad (A,B) il problema
%  -MU*D(DU/DX)/DX + ETA*DU/DX + SIGMA*U = BVPFUN
%  sull'intervallo (A,B) con condizioni al bordo
%  U(A)=UA e U(B)=UB. BVPFUN puo' essere un function
%  handle o una user defined function.
%  XH e UH contengono, rispettivamente, i nodi
%  e la soluzione numerica, inclusi i valori al bordo.
h = (b-a)/(N+1);
xh = (linspace(a,b,N+2))';
hm = mu/h^2;
hd = eta/(2*h);
e =ones(N,1);
Afd = spdiags([-(hm+hd)*e (2*hm+sigma)*e (-hm+hd)*e],...
    -1:1, N, N);
xi = xh(2:end-1);
f =bvpfun(xi,varargin{:});
f(1) =  f(1)+ua*(hm+hd);
f(end) = f(end)+ub*(hm-hd);
uh = Afd\f;
uh=[ua; uh; ub];
