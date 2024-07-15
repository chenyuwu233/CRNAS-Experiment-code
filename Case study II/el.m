

function ret = el(theta,DATA,Time,Init)
    p = theta(1);
    r2 = theta(2);
    m1 = theta(3);
    Init = [p*Init,(1-p)*Init];

    f_el = @(t,y) [y(1)*(1-(y(1)+r2*y(2))/(m1));r2*y(2)*(1-(y(1)+r2*y(2))/(m1))];
    
    [~,DATA_est] = ode45(f_el,Time,Init);

    if size(DATA_est) ~= size(DATA)
        size(DATA)
        pause
    end
    ret = sum((DATA-DATA_est).^2,'all');


end