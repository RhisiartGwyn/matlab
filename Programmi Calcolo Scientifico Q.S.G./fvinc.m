function [f]=fvinc(t,y)
% FVINC Function per  l'esempio del pendolo sferico
[n,m]=size(y); f=zeros(n,m);
phix=2*y(1); phiy=2*y(2); phiz=2*y(3);
H=2*eye(3); mass=1;
xp=zeros(3,1); xp(1:3)=y(4:6);
F=[0;0;-mass*9.8];
G=[phix;phiy;phiz];
lambda=(mass*xp'*H*xp+F'*G)/(G'*G);
f(1:3)=y(4:6);
for k=1:3;
  f(k+3)=(F(k)-lambda*G(k))/mass;
end
