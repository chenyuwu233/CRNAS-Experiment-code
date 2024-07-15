% This function generate the data for the whole time vector so the return
% will be 1xNT vector


function ret = Sto_samplepath(initial,theta,d,Time,pi)
    ret = zeros(1,length(Time));
    X   = round(initial.*pi);
    B_D   = zeros(2,size(theta,2));
    for i = 1:size(theta,2)
        B_D(1,i) = theta(1,i);
        B_D(2,i) = theta(2,i) - log(theta(3,i) + (1 - theta(3,i))/(1 + (d/theta(4,i))^(theta(5,i))));
    end
    a      = B_D(1,:) + B_D(2,:);
    t      = Time(1);
    ret(1) = initial;
    nt     = 2;
    while t < Time(length(Time)) && sum(X) >0
       rate   = X * a';
       rate_v = X.*a;
       t      = t - (1/rate) * log(rand);
       if t > Time(nt)
           ret(nt) = sum(X);
           nt      = nt+1;
           if nt > length(Time)
               break
           end
       end
       if sum(X) > 3e6
           break
       end
       
       event  = rand;
       for i = 1:size(rate_v,2)
          block = sum(rate_v(1:i))/rate;
          if event <= block
              event_bd = rand;
              if event_bd <= B_D(1,i)/(B_D(1,i) + B_D(2,i))
                  X(i)      = X(i) + 1;
              else
                  X(i)      = X(i) - 1;
              end
              break 
          end
       end
    end
end