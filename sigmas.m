function [condition,test,x] = sigmas(T,mu,s) %wrapper program for solving chiral field 
cla;
numMesh = 1000;
ui = 1e-2;
% uf = 1-(1e-2);
uf = 9.8e-1;
u = linspace(ui,uf,numMesh);

mu0    = 430;
mu1    = 830;
mu2    = 176;
% v3     = -3;
v3=0;   %two flavor
v4     = 8;
mq     = 0*9.75;
%kappa  = 1; %defined in bhsolve


[zh,Q] = bhsolve(T,mu);

f=1-(1 + Q^2)*(u.^4)+(Q^2)*(u.^6);
df=gradient(f,u);

[chi]=chisolve(u,zh,Q,mu0,mu1,mu2,v3,v4,mq,s);

ch = transpose(chi);
x=ch(1,:); %chiral field
dx=ch(2,:); %derivative of chiral field

test = (-(u.^2).*df).*dx./f - ((3*x-3*v3*(x.^2))-4*v4*(x.^3))./f;

% plot(u,chi(:,1));
% plot(u,test);

% compare dilaton derivatives
% phi=-(mu1^2*zh^2)*(u.^2)+(mu1^2+mu0^2)*(zh^2)*(u.^2).*(1-exp(-(((mu2*zh)^2)*(u.^2))));
% dph=gradient(phi,u);
% dphn= -2*zh^2*mu1^2*u+2*zh^2*(mu0^2+mu1^2).*u.*(1-exp(-(((mu2*zh)^2)*(u.^2))))+2*mu2^2*zh^4*u.^3*(mu0^2+mu1^2).*exp(-(((mu2*zh)^2)*(u.^2)));
% plot(u,dph);
% plot(u,dphn);


condition=test(end-2);


end
