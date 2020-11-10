syms f(x) 
f(x) = sin(x);
df = diff(f);
format long
x0 = 2;
N = 30;
x = x0;
k = 0;
tol = 1e-15;
err = 1;

while k < N && err > tol
    
    x1 = x-(f(x)/df(x));
    double (x1)
    err = abs(double(x) - double(x1))
    x = x1;
    k = k +1;
    
end

err < tol
