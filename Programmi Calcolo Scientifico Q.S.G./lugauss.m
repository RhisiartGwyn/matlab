function A=lugauss(A)
%LUGAUSS   Fattorizzazione LU senza pivoting
%   A = LUGAUSS(A) calcola la fattorizzazione
%   LU della matrice A, memorizzando nella
%   parte triangolare inferiore stretta di A la
%   matrice L (gli elementi diagonali di L sono tutti
%   uguali a 1) ed in quella superiore il fattore U
[n,m]=size(A);
if n ~= m;
 error('A non e'' una matrice quadrata'); else
 for k = 1:n-1
   for i = k+1:n
     A(i,k) = A(i,k)/A(k,k);
     if A(k,k) == 0;
     error('Un elemento pivot si e'' annullato');
     end
     j = [k+1:n]; A(i,j) = A(i,j) - A(i,k)*A(k,j);
   end
 end
end
