function [zh, Q] = bhsolve(T,mu)

kappa=1;

if (mu > 0)
    lb = [0,0]; % lower bound constraint
    rng default % reproducible initial point
    x0 = zeros(2,1); %zh and q start at 0
    f = @(x)fminconstr(x,T,mu);
    opts = optimoptions(@fmincon,'Algorithm','interior-point','Display','off');
    x = fmincon(@(x)0,x0,[],[],[],[],lb,[],f,opts);
    zh=x(1);
    q=x(2);
else
    zh=1/(pi*T);
    q=0;
end
    
    Q=q*zh^3;
    %temp = (1-(q^2)*(zh^6)/2)/(pi*zh); %check temperature
    %bcp =  kappa*q*zh^2; %check baryon chemical potential
   

%==========================================================================    
% optimization functions
    
    
    function F = fbh(x,T,mu) %define the equations for finding zh and q
        F(1) = T*(pi*x(1)) - (1-(x(2)^2)*(x(1)^6)/2);
        F(2) = mu - kappa*x(2)*x(1)^2;
    end
    
    function [c,ceq] = fminconstr(x,T,mu)
        c = []; % no nonlinear inequality
        ceq = fbh(x,T,mu); % the fsolve objective is fmincon constraints
    end
    

end

