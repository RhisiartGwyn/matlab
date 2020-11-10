syms f(x) 
f(x) = x^2-2;
df = diff(f);
format long
x0 = 1;
N=5;
x= x0;

for k = 1:N
    
    x = x-(f(x)/df(x));
    double (x)
end
