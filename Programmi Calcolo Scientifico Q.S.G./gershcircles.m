function gershcircles(A)
%GERSHCIRCLES disegna i cerchi di Gershgorin
%   GERSHCIRCLES(A) disegna i cerchi di
%   Gershgorin associati alla matrice  A ed
%   alla sua trasposta.
n = size(A);
if n(1) ~= n(2)
    error('Solo matrici quadrate');
else
    n = n(1);
    circler = zeros(n,201);
    circlec = circler;
end
center = diag(A);
radiic = sum(abs(A-diag(center)));
radiir = sum(abs(A'-diag(center)));
one = ones(1,201);
cosisin = exp(i*[0:pi/100:2*pi]);
figure(1); title('Cerchi riga','Fontsize',20);
set(gca,'Fontsize',20)
xlabel('Re'); ylabel('Im');
figure(2); title('Cerchi colonna','Fontsize',20);
set(gca,'Fontsize',20)
xlabel('Re'); ylabel('Im');
color=[0.8,1,1];
for k = 1:n
   circlec(k,:) = center(k)*one + radiic(k)*cosisin;
   circler(k,:) = center(k)*one + radiir(k)*cosisin;
   figure(1);
   patch(real(circler(k,:)),imag(circler(k,:)),...
   color(1,:));
   hold on
   plot(real(circler(k,:)),imag(circler(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
   figure(2);
   patch(real(circlec(k,:)),imag(circlec(k,:)),...
   color(1,:));
   hold on
   plot(real(circlec(k,:)),imag(circlec(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
end
for k = 1:n
   figure(1);
   plot(real(circler(k,:)),imag(circler(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
   figure(2);
   plot(real(circlec(k,:)),imag(circlec(k,:)),'k-',...
   real(center(k)),imag(center(k)),'kx');
end
