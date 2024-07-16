function ret = mix_logistic_nl_obj(theta, DATA,cv ,inputs)
    init = inputs.init;

    Theta = reshape(theta,3,[]);
    num_sub = size(Theta,2);
    est = 0;
    for i = 1:num_sub
        est = est + init*Theta(1,i)*logistic(cv,1,Theta(2,i),Theta(3,i));
    end
    ret = (est-DATA)^2;
end