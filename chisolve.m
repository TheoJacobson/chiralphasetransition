function [y] = chisolve(u,zh,Q,mu0,mu1,mu2,v3,v4,mq,s)

sig=s^3;

initialC = [mq*sqrt(3)*zh*u(1)/(2*pi) + sig*2*pi*(zh^3)*u(1)^3/(sqrt(3))
              mq*sqrt(3)*zh/(2*pi) + sig*6*pi*(zh^3)*u(1)^2/(sqrt(3))];
 
           options = odeset('RelTol',1e-10,'AbsTol',1e-9);
       
    try
           [~,y]=ode45(@(u,y) dydx(u,y),u,initialC,options);

    catch
            disp(s);
    end
    function div = dydx(u,y)
        
     f=1-(1+Q^2)*(u^4)+(Q^2)*(u^6);
     dphi= -2*zh^2*mu1^2*u+2*zh^2*u*(mu0^2+mu1^2)*(1-exp(-(((mu2*zh)^2)*(u^2)))) ...
              +2*mu2^2*zh^4*u^3*(mu0^2+mu1^2)*exp(-(((mu2*zh)^2)*(u^2)));
     
     eq1 = y(2); 
     eq2 = (4-f+u*f*dphi)*y(2)/(u*f)-(3*y(1)-3*v3*y(1)^2-4*v4*y(1)^3)/((u^2)*f);
         
     div=[eq1; eq2];
     
    end
 

end
