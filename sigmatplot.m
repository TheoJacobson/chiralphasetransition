function [sigmavst] = sigmatplot(mu,Tmin,Tmax,numT) %solves for sigma for various temperatures

T=Tmin;
dT = (Tmax-Tmin)/numT;
sigmavst = zeros(2,numT);
i=1;

while T < Tmax
    sigmavst(1,i) = T;
    sigmavst(2,i) = (sigmasolve(T,mu)*0.001)^3; % MeV^3
    i = i+1;
    T = T+dT;
    disp(T);
end

cla;
plot(sigmavst(1,:),sigmavst(2,:));

end