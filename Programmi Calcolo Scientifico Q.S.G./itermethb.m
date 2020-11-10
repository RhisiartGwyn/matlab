function [x,iter,e]= itermethb(B,g,x0,tol,kmax)
% ITERMETHB metodo iterativo classico per sist. lineari
% [X,ITER,E]= ITERMETHB(B,G,X0,TOL,KMAX)
% Calcola la soluzione del metodo iterativo
% x^{k+1}=B x^{k}+G, con x^{0} assegnato, per k=0,1....
% B e' la matrice di iterazione, il vettore iniziale
% e' memorizzato in X0. TOL e KMAX sono tolleranza e
% numero massimo di iterazioni per il test d'arresto
% sull'incremento. La soluzione e' memorizzata in X,
% ITER e' il numero di iterazioni richieste per
% giungere a convergenza.
% E e' il vettore degli errori ||x^{k+1}-x^{k}||.
k=0;
err=tol+1;
e=[];
while err> tol && k< kmax
    x=B*x0+g;
    err=norm(x-x0);
    e=[e;err];
    x0=x;
    k=k+1;
end
iter=k;
