function [y,b] = horner(a,z)
%HORNER Metodo di Horner
%   Y=HORNER(A,Z) calcola
%   Y = A(1)*Z^N + A(2)*Z^(N-1) + ... + A(N)*Z + A(N+1)
%   con il metodo di divisione sintetica di Horner.
n = length(a)-1; b = zeros(n+1,1);
b(1) = a(1);
for j=2:n+1
   b(j) = a(j)+b(j-1)*z;
end
y = b(n+1); b = b(1:end-1);
