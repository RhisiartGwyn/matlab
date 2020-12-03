v = input('Vettore: ');
p =1; 
m = v(1);
k = 2;

for k = 2:length(v)
    if v(k) >= m
        m = v(k);
        p = k;
    end
end
