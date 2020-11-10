function [x,err,k]=penalty(fun,grad_fun,h,grad_h,...
g,grad_g,x0,tol,kmax,kmaxd,meth,varargin)
% PENALTY Ottimizzazione vincolata con penalizzazione
%  [X,ERR,K]=PENALTY(FUN,GRAD_FUN,H,GRAD_H,...
%  G,GRAD_G,X0,TOL,KMAX,KMAXD,METH)
%  Approssima un punto di minimo della funzione FUN
%  soggetto ai vincoli di uguaglianza H=0 e di disu-
%  guaglianza G>=0, con un metodo di penalizzazione,
%  partendo da un punto X0, fino al
%  raggiungimento della tolleranza TOL o di KMAX
%  iterazioni. GRAD_FUN, GRAD_H, e GRAD_G contengono
%  i gradienti di FUN, H e G, rispettivamente. I
%  vincoli ed i rispettivi gradienti sono da
%  inizializzare con [], qualora non siano presenti.
%  Per risolvere il problema di minimo non vincolato
%  si richiama il programma FMINSEARCH (se METH=0)
%  o DESCENT (se METH>0). Quando METH>0, KMAXD e
%  METH contengono rispettivamente il numero
%  massimo di iterazioni e la scelta del metodo
%  di discesa per il programma DESCENT. Se METH=1
%  (o METH=2) si richiede l'espressione dell'Hes-
%  siana (o una sua approssimazione al passo 0) come
%  ultimo parametro in input.
xk=x0(:); alpha0=1;
if meth==1, hess=varargin{1};
elseif meth==2, hess=varargin{1};
else  hess=[]; end
if ~isempty(h), [nh,mh]=size(h(xk)); end
if ~isempty(g), [ng,mg]=size(g(xk)); else, ng=[]; end
err=tol+1; k=0;
alphak=alpha0; alphak2=alphak/2; told=.1;
while err>tol && k< kmax
  P=@(x)Pf(x,fun,g,h,alphak2,ng);
  grad_P=@(x)grad_Pf(x,grad_fun,h,g,grad_h,...
                     grad_g,alphak,ng);
  if meth==0
     options=optimset('TolX',told);
     [x,err,kd]=fminsearch(P,xk,options);
     err=norm(x-xk);
  else
     [x,err,kd]=descent(P,grad_P,xk,told,kmaxd,...
                        meth,hess);
     err=norm(grad_P(x));
  end
  if kd<kmaxd
    alphak=alphak*10; alphak2=alphak/2;
  else
    alphak=alphak*1.5; alphak2=alphak/2;
  end
  k=k+1; xk=x; told=max([tol,told/10]);
  end
end % end of the function penalty
function y=Pf(x,fun,g,h,alphak2,ng)
y=fun(x);
if ~isempty(h), y=y+alphak2*sum((h(x)).^2); end
if ~isempty(g), G=g(x);
for j=1:ng, y=y+alphak2*max([-G(j),0])^2; end
end
end % end of function Pf
function y=grad_Pf(x,grad_fun,h,g,...
                  grad_h,grad_g,alphak,ng)
y=grad_fun(x);
if ~isempty(h), y=y+alphak*grad_h(x)*h(x); end
if ~isempty(g), G=g(x); Gg=grad_g(x);
for j=1:ng
if G(j)<0, y=y+alphak*Gg(:,j)*G(j); end
end, end
end % end of function grad_Pf
