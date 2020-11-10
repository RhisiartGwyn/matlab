%input function using anonymous function syntax 
f = input('Introduci la funzione f = ')
%input interval boundaries
a = input('Introduci l''estremo inferiore dell''intervallo a = ')
b = input('Introduci l''estremo superiore dell'' intervallo b = ')  %#ok<*NOPTS>
%input error (or tolerance)
prec = input ('Precisione = ')
%compute number of iterations needed, rounded towards plus infinity
k = ceil(log2(b-a)-log2(prec));

%compute midpoint
xpm = (a + b)/2;
%evaluate function
ypm = f(xpm);
%print results for k = 1
fprintf ('k = 1 \t a = %16.15f \t b = %16.15f \t ', a, b) 
fprintf ('x_k = %16.15f \t f(x_k) = %16.15f\n', xpm, ypm)

for i = 2:k                                 %start algorithm cycle  
    
    if f(a) * ypm < 0                       %if true, root in (a, xpm) interval
        b = xpm;

    else                                    %if not, root in (xpm, b) interval
        a = xpm;

    end
    
    xpm = (a + b)/2;                        %midpoint of previously found interval
    ypm = f(xpm);                           %evaluate function at new midpoint
    
    %print results for each cycle k
    fprintf ('k = %d \t a = %16.15f \t b = %16.15f \t ', i, a, b)
    fprintf ('x_k = %16.15f \t f(x_k) = %16.15f\n', xpm, ypm)
    
end

%find zero of the function using matlab's own tool (16 digits precision), 
%taking bisection method's result as starting point, and makes a comparison
fzero (f, xpm);
fprintf ('\nErrore: %16.15f\n', abs (fzero (f, xpm) - xpm))

fplot (f)
hold on
axis ([a-0.1 b+0.1 -0.1 0.1])
xline (xpm)
yline (0)
