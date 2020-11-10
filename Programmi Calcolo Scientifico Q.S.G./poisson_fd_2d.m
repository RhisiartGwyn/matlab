function [xh,yh,uh,error]=poisson_fd_2d(a,b,c,d,...
                         nx,ny,fun,gd,uex,varargin)
%POISSON_FD_2D approssimazione del problema di Poisson
%   in due dimensioni
%  [XH,YH,UH]=POISSON_FD_2D(A,B,C,D,NX,NY,FUN,GD)
%  risolve con lo schema alle differenze finite
%  a 5 punti il problema -LAPL(U) = FUN in un
%  rettangolo (A,B)X(C,D)  con condizioni al bordo
%  di Dirichlet U(X,Y)=GD(X,Y) per ogni (X,Y)
%  sul bordo del rettangolo.
%  [XH,YH,UH,ERROR]=POISSONFD(A,B,C,D,NX,NY,FUN,...
%  GD,UEX) calcola anche l'errore sulla soluzione
%  esatta UEX.
%  FUN,GD e UEX possono function handle o user
%  defined functions.
%  [XH,YH,UH,ERROR]=POISSON_FD_2D(A,B,C,D,NX,NY,FUN,...
%  GD,UEX,P1,P2, ...) passa i parametri opzionali
%  P1,P2,... alle funzioni FUN,GD,UEX.
if nargin == 8
    uex = @(x,y) 0*x+0*y;
end
nx1 = nx+2; ny1=ny+2; dim = nx1*ny1;
hx = (b-a)/(nx+1); hy = (d-c)/(ny+1);
    hx2 = hx^2;      hy2 = hy^2;
aii = 2/hx2+2/hy2; aix = -1/hx2;  aiy = -1/hy2;
A = speye(dim,dim);  f = zeros(dim,1);
y = c;
for m = 2:ny+1
 x = a; y = y + hy;
 for n = 2:nx+1
   i = n+(m-1)*nx1; x = x + hx;
   f(i) = fun(x,y,varargin{:});
   A(i,i) = aii; A(i,i-1) = aix; A(i,i+1) = aix;
   A(i,i+nx1) = aiy;   A(i,i-nx1) = aiy;
 end
end
ub = zeros(dim,1); xh = [a:hx:b]'; yh = [c:hy:d];
ub(1:nx1) = gd(xh,c,varargin{:});
ub(dim-nx-1:dim) = gd(xh,d,varargin{:});
ub(1:nx1:dim-nx-1) = gd(a,yh,varargin{:});
ub(nx1:nx1:dim) = gd(b,yh,varargin{:});
f = f - A*ub;
nbound = [[1:nx1],[dim-nx-1:dim],[1:nx1:dim-nx-1],...
    [nx1:nx1:dim]];
ninternal = setdiff([1:dim],nbound);
A = A(ninternal,ninternal);
f = f(ninternal);
utemp = A\ f;
u = ub; u (ninternal) = utemp;
k = 1; y = c;
for j = 1:ny1
    x = a;
    for i = 1:nx1
        uh(j,i) = u(k);         k = k + 1;
        ue(j,i) = uex(x,y,varargin{:});
        x = x + hx;
    end
    y = y + hy;
end
if nargout == 4 && nargin >= 9
    error = max(max(abs(uh-ue)))/max(max(abs(ue)));
elseif nargout == 4 && nargin ==8
    warning('Soluzione esatta non disponibile');
    error = [ ];
end
end
