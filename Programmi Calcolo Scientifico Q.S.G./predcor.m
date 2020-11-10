function [t,u]=predcor(odefun,tspan,y0,Nh,...
               predictor,corrector)
%PREDCOR  Risolve equazioni differenziali
%  con il metodo predictor-corrector
%  [T,Y]=PREDCOR(ODEFUN,TSPAN,Y0,NH,@PRED,@CORR) con
%  TSPAN=[T0 TF] integra il sistema di equazioni
%  differenziali y' = f(t,y) dal tempo T0 a TF con
%  condizioni iniziali Y0 utilizzando un generico
%  metodo predictor-corrector su una griglia
%  equispaziata di NH intervalli.
%  La funzione ODEFUN(T,Y) deve ritornare un vettore
%  contenente f(t,y), della stessa dimensione di y.
%  Ogni riga del vettore soluzione Y corrisponde ad
%  un istante temporale del vettore colonna T.
%  Le function PRED e CORR caratterizzano il tipo
%  di metodo predictor-corrector scelto.
%  Possono essere scelte, ad esempio, tra le
%  function proposte:
%  feonestep, beonestep, cnonestep.
h=(tspan(2)-tspan(1))/Nh;
t=linspace(tspan(1),tspan(2),Nh+1)';
y=y0(:); % genera sempre un vettore colonna
d=length(y0);
u=zeros(Nh+1,d);
u(1,:)=y0.'; % trasposta anche di variabili complesse
for n=1:Nh
    wn=u(n,:).';
    fn = odefun(t(n),wn);
    upre = predictor(t(n),wn,h,fn);
    w = corrector(t(n+1),wn,upre,h,odefun,fn);
    u(n+1,:)=w.';
end
end
