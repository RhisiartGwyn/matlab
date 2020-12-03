clear all

n = input('Inserisci un vettore di 4 numeri interi e positivi: ');

if ~all(n > 0) || length(n) ~= 4 || ~all(n == round(n))
    error('Baka!')
end

x = linspace (-1, 1, 100);
x1 = x(x <= 0);
x2 = x(x > 0);

figure

for i = 1:4
    f1 = x1.^n(i) + 1;
    f2 = (-1).^n(i)*abs(x2).^n(i);
    subplot(2, 2, i)
    plot(x1, f1, x2, f2)
    title(['grafico per n = ', int2str(n(i))])
end
