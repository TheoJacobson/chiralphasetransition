function [suscmax,out] = sigmatplot(mq,mu,Tmin,Tmax,dT) %solves for sigma for various temperatures

axis([0.06 Tmax*0.001 0 0.13]);
T=Tmin;
numT = 500; % big number
s0=0;
dTmin = 0.5; % minimum T step
sigmavst = zeros(2,numT);
prerun = zeros(2,(Tmax-Tmin)/dT);
i=1;
check=0;

%========================================== check for max derivative
while T < Tmax
    
    prerun(1,i) = T;
    [s] = sigmasolve(mq,T,mu);
    prerun(2,i) = s;
    i = i+1;
    T = T+dT;
end
    dcrit = min(gradient(prerun(2,:),prerun(1,:)));
    
%========================================== start run-through
T=Tmin;
i=1;

    while (dT >= dTmin)
        if(check==1)
        Tmax=Tc+3*dT;
        T=Tc-3*dT;
        dT=0.49*dT;
        end
        check=0;
    
while T < Tmax
    sigmavst(1,i) = 0.001*T;
    [s] = sigmasolve(mq,T,mu);
    
    if(check==0)
        if((s-s0)/dT <= dcrit)
         Tc=T;
         check=1;
         disp(Tc);
        end
    end
    
    s0=s;
    sigmavst(2,i) = (s*0.001)^3; % MeV^3
    i = i+1;
    T = T+dT;
end
    end

%=============================================== plots / susceptibility
%cla;
out = sortrows(sigmavst',1)'; % sorts by temp
out = out(:,any(out)); % removes extra columns
susc = -gradient(out(2,:),out(1,:));
suscmax = max(susc);

% xlabel('T (GeV)');
% ylabel('\sigma (GeV^3)');
% plot(out(1,:),out(2,:));
% scatter(out(1,:),out(2,:),'.');

% plot(out(1,:),susc);


end

