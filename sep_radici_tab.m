% sep_radici_tab2
% script per la separazione per tabulazione delle radici dell’equazione non
% lineare f(x)=0 contenute nell’intervallo I = [a,b]
% Lo script richiede in input:
% f = espressione della funzione di cui separare gli zeri in forma di anonymous function
% a = estremo sinistro dell’intervallo (dominio) in cui si vogliono separare gli zeri
% b = estremo destro dell’intervallo (dominio) in cui si vogliono separare gli zeri
% np = numero di punti equidistanti nell’intervallo [a,b] da usare per la tabulazione
%
% Lo script stampa e salva nella variabile intervallo gli estremi degli intervalli di
% separazione.
% La variabile intervallo è una matrice la cui prima riga contiene gli estremi inferiori degli
% intervalli di separazione, mentre la seconda riga contiene i corrispondenti estremi superiori

clear all
f = input('Funzione: ');

if ~isa(f, 'function_handle')
    error('Baka!')
end

a = input('Estremo inferiore: ');
b = input('Estermo superiore: ');

if a >= b
    error('Baka!')
end

np = input('Punti di suddivisione dell''intervallo: ');

if np <= 0
    error('Baka!')
end

x = linspace (a, b, np);
y = f(x);
pos = find(abs(diff(sign(y))));

if isempty(pos)
    error('Baka!')
end

intervallo(1,:) = x(pos);
intervallo(2,:) = x(pos+1);

fprintf('estremo inferiore \t estremo superiore \n')
fprintf('%16.15f \t %16.15f\n', intervallo)

plot(x, y)
yline(0)