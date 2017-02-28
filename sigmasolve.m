function [egnvalues] = sigmasolve(mq,T,mu) %solves for sigma 

numvalues=1;
egnvalues=zeros(1,numvalues);    % list of eigenvalues to be filled
reset=100;                       % starting step size
minstep=0.5;                   % smallest step size                   
condition=0;
count=0;
s=500;  % initial energy level
% smax=0;
% safety=0;

while count < numvalues   % while we fill up the eigenvalue list
    
smax=0;
safety=0;
    
   step=reset;
  while step>=minstep    
    oldcondition = condition;
    [condition] = sigmas(mq,T,mu,s);
    
    if (oldcondition*condition) < 0
        s = s + step;
        step = step*0.1;
        condition = oldcondition;
    end
    s = s - step;
    if s < smax
        safety=1;
       break;
    end
  end
    count = count+1;
    if safety==1
      egnvalues(count)=0;
    else
        egnvalues(count) = s;
    end
    
%     if safety==1
%       disp('SAFETY');    
%        break
%     end
    
    step = reset; 
%     [~,chi,test] = sigmas(T,mu,s);
%     plot(chi);
%     plot(test);
    condition = 0;
    s = s-step;
    

end

end
