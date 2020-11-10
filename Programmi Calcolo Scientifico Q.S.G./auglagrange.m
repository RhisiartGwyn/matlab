function [x,err,k]=auglagrange(fun,grad_fun,h,grad_h,...
    x0,lambda0,tol,kmax,kmaxd,meth,varargin)
% AUGLAGRANGE Ottimizz. vincolata con Lagr. aumentata
%  [X,ERR,K]=AUGLAGRANGE(FUN,GRAD_FUN,H,GRAD_H,...
%  X0,LAMBDA0,TOL,KMAX,KMAXD,METH)
%  Approssima un punto di minimo della funzione FUN
%  soggetto ai vincoli di uguaglianza H con il metodo
%  della Lagrangiana aumentata, partendo da un punto
%  X0, fino al raggiungimento della tolleranza TOL
%  o di KMAX iterazioni. GRAD_FUN e GRAD_H contengono
%  i gradienti di FUN e H, rispettivamente.
%  Per risolvere il problema di minimo non vincolato
%  si richiama il programma FMINSEARCH (se METH=0)
%  o DESCENT (se METH>0). Quando METH>0, KMAXD e
%  METH contengono rispettivamente il numero
%  massimo di iterazioni e la scelta del metodo
%  di discesa per il programma DESCENT. Se METH=1
%  (o METH=2) si richiede l'espressione dell'Hessiana
%  (o una sua approssimazione al passo 0) come
%  ultimo parametro in input.
alpha0=1;
if meth==1, hess=varargin{1};
elseif meth==2, hess=varargin{1};
else, hess=[]; end
err=tol+1; k=0; xk=x0(:); lambdak=lambda0(:);
if ~isempty(h), [nh,mh]=size(h(xk)); end
alphak=alpha0; alphak2=alphak/2; told=0.1;
while err>tol && k< kmax
  L=@(x)Lf(x,fun,lambdak,alphak2,h);
  grad_L=@(x)grad_Lf(x,grad_fun,lambdak,alphak,...
                     h,grad_h);
  if meth==0
    options=optimset('TolX',told);
    [x,err,kd]=fminsearch(L,xk,options);
    err=norm(x-xk);
  else
    [x,err,kd]=descent(L,grad_L,xk,told,kmaxd,meth,...
                       hess);
    err=norm(grad_L(x));
  end
  lambdak=lambdak-alphak*h(x);
  if kd<kmaxd
    alphak=alphak*10; alphak2=alphak/2;
  else
    alphak=alphak*1.5; alphak2=alphak/2;
  end
  k=k+1; xk=x; told=max([tol,told/10]);
  end
end % end auglagrange

function y=Lf(x,fun,lambdak,alphak2,h)
y=fun(x);
if ~isempty(h)
  y=y-sum(lambdak'*h(x))+alphak2*sum((h(x)).^2);
end
end % end function Lf

function y=grad_Lf(x,grad_fun,lambdak,alphak,h,grad_h)
y=grad_fun(x);
if ~isempty(h)
   y=y+grad_h(x)*(alphak*h(x)-lambdak);
end
end % end function grad_Lf
