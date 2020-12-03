clear all

m = input('Inserisci un vettore di 6 numeri interi e positivi: ');

if ~all(m > 0) || length(m) ~= 6 || ~all(m == round(m))
    error('Baka!')
end

n = input('Inserisci il numero di punti in cui suddividere l''intervallo: ');

if ~all(n > 0)
    error('Baka!')
end

figure

for i = 1:6   
    x = linspace (-m(i), m(i), n);
    x1 = x(x <= 0);
    x2 = x(x > 0);
    
    axis tight
    f1 = (m(i) - (x1.^2)/m(i)).^m(i);
    f2 = ((x2.^2)/m(i) + m(i)).^m(i);
    subplot(3, 2, i)
    plot(x1, f1, x2, f2)
    title(['grafico per n = ', int2str(m(i))])
end
