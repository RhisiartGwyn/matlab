% Lo script esercizio_lab2 legge in input il numero intero x0, grafica
% i primi 500 elementi della successione x(n) e ne stampa il massimo
% e il minimo.
% La successione x(n) è definita per n>=1 come segue
% x(n)= 2 x(n-1) + 2 se n è dispari
% x(n)= x(n-1)/2 + 1 se n è pari

x0 = input('Inserisci il punto iniziale (numero intero): ');

if ((length(x0)>1) || (x0~=round(x0)))
    error('Il punto iniziale deve essere un numero intero!!!')
end

x(1) = x0;

for n=2:500
    if rem(x(n-1),2)==0
        x(n) = x(n-1)/2 + 1;
    else
        x(n) = 2*x(n-1) + 2;
    end
end

figure
plot(x)
title(['Successione con punto iniziale ', int2str(x0)])

M = max(x);
m = min(x);

fprintf('Il massimo e il minimo della successione sono %10d \t %10d\n', M, m)