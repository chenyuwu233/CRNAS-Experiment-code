function ret = test_obj(theta, DATA,cv ,inputs)
    inputs.cmd = 0;
    ret = mix_logistic_nl_full(theta,DATA,cv,inputs);
end