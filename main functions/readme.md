In this repository, we include main functions(CRNAS.m, cub_new_sub.m, mynonconvexqp.m) for CRNAS implementation and several example test code (main_mix_logistic_test.m).
Those functions require several types of inputs:
* function input: for this type of input we will specify the function format
* scalar input
* vector input
* structure input

* [x, val,jj] = CRNAS(func_all,func_obj,func_inputs,x,DATA,Design,matrix_A,vec_b,lb,ub,max_iter,opt_tol,step_tol): main structure for CRNAS algorithm.
    * Inputs: 
        * func_all (function input): function that return objective value, gradient, and Hessian of the objective function.
            Example: function [obj,grad,Hess] = func_all(theta, DATA, Design, inputs).
            theta: targeted parameters.
            DATA: dataset of interest.
            Design: corresponding designed variables, e.g. Time when estimating the cell population dynamics.
            inputs: fixed inputs related to the objective computation(usualy in a structure format).
        * func_obj (function input): function that only return objective value.
            Example: function obj = func_obj(theta,DATA,Design,inputs).
        * func_inputs (structure input): structure of all information necessary to compute the objective function.
        * x (vector inputs): starting point of the optimization process.
        * DATA (vector inputs): dataset of interest.
        * Design (vector inputs): corresponding designed variables, should have a same dimension as DATA.
        * matrix_A (vector inputs): matrix A in equality constraints Ax = b.
        * vec_b (vector inputs): vector b in equality constraints Ax = b.
        * lb (vector inputs): lower bound for variable x.
        * ub (vector inputs): upper bound for variable x.
        * max_iter (scalar input): maximum number of iteration.
        * opt_tol (scalar input): tolerance for the first-order optimality condition.
        * step_tol (scalar input): tolerance for the step size.
    * Outputs:
        * x: optimal solution
        * val: optimal objective value
        * jj: exiting outer iteration number 
