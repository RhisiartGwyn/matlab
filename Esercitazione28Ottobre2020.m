%Esercizio 1
%Porre a zero gli elementi negativi di x e assegnare il risultato alla variabile y:
x = [7 6 1 2 0 -1 4 3 -2 0];
y = x;
y(y < 0) = 0;

%Generare il vettore t contenente gli elementi di x più grandi di 3:
t = x(x > 3);

%Aggiungere 3 agli elementi di indice pari di x e assegnare il risultato alla variabile z:
z = x;
z(2:2:end) = z(2:2:end) + 3;

%Detto m il valore della media degli elementi di z:
m = mean (z);

%Assegnare il valore 0 agli elementi di z inferiori a m:
z(z < m) = 0;

%Agli elementi di z superiori alla media, assegnare il valore della differenza degli stessi con m:
z(z > m) = z(z > m) - m;

%Aggiungere 3 agli elementi pari di x e assegnare il risultato alla variabile w:
w = x + 3 * (rem (x, 2) == 0);

%Detto m2 il valore della media degli elementi di w, assegnare il valore 0 agli elementi di w inferiori a m2:
m2 = mean (w);
w(w < m2) = 0;

%Assegnare alla variabile p il valore del prodotto degli elementi di t:
p = prod (t);

%Assegnare alla variabile A la matrice le cui righe sono, nellassegnare alla variabile A la matrice le cui righe sono, nell'ordine, i vettori y, z, w:
A = [y; z; w];


%Assegnare alla variabile M il valore della media degli elementi di A:
M = mean (A, 'all');


%Assegnare alla variabile P il valore del prodotto degli elementi di A:
P = prod (A, 'all');






%Esercizio 2
%Scrivere i comandi Matlab necessari per:
% 1. separare la radice dell’equazione non lineare associata al problema usando il metodo della tabulazione.
% 2. stampare su un’unica riga il messaggio l'estremo inferiore dell’intervallo di separazione e seguito dal valore determinato.
% 3. stampare su un’unica riga il messaggio l'estremo superiore dell’intervallo di separazione e seguito dal valore determinato.
% 4. sulla stessa finestra grafica, tracciare il grafico della funzione di cui si cerca lo zero e dell asse delle ascisse, usando tratti di linea e colori distinti. Il grafico deve contenere anche gli estremi dell intervallo di separazione determinato --- si usi un puntino nero come marcatore di punto di dimensione 5. Il grafico deve essere completo    di etichette per gli assi, titolo e legenda.            
% 5. salvare la figura nel file figuraesercizio2.png5.

u = 2510; g = 9.81; M_0 = 2.8e+06; m = 13.3e+3; a = 335; prec = 1e-3;
v = @(t)(u .* log(M_0 ./ (M_0 - m .* t)) - g .* t);
t = (0:prec:100);
y = v(t) - a;
pos = find(diff(sign(y)));
disp (['L''estremo inferiore dell''intervallo di separazione è ' num2str(t(pos))])
disp (['L''estremo superiore dell''intervallo di separazione è ' num2str(t(pos+1))])
axis([t(pos)-prec t(pos+1)+prec -prec prec])
hold on
plot (t, y, 'r--', t, zeros(1, length(t)), 'k')
plot (t(pos), 0, 'k.', t(pos+1), 0, 'k.', "LineWidth", 5)
title ('Titolo')
xlabel ('Asse x')
ylabel ('Asse y')
legend ('v(t) - a')
saveas (gcf, 'figuraesercizio2', 'png')
