function [suscmq] = mqdependence(mqmin,mqmax,dmq)

num = floor((mqmax-mqmin)/dmq);
suscmq = zeros(2,num);

mu=0;

mq=mqmin;
i=1;

while mq < mqmax
   suscmq(1,i) = mq;
   [susc,~] = sigmatplot(mq,mu,120,220,10);
   suscmq(2,i) = susc;
    i=i+1;
    mq=mq+dmq;
   disp(mq);
end

plot(suscmq(1,:),suscmq(2,:));


end