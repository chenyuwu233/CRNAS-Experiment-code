function [x, val,jj]=CRNAS(func_all,func_obj,x,matrix_A,vec_b,lb,ub,max_iter,opt_tol,step_tol)
    step_size=6e-1;
    thro=step_size^2;
    if any(lb>ub)
        error('upper bound is lower than lower bound')
    end
    if any(x<lb) || any(x>ub)
        error('x is out of bound')
    end
    if any((x==lb)&(x~=ub)) || any((x~=lb)&(x==ub))
        error('x is not an interior point')
    end
    if norm(matrix_A*x-vec_b)>1e-6
        keyboard
        error('x is not feasible')
    end 
    set_equal=(lb==ub);
    lb(set_equal)=-inf;
    ub(set_equal)=inf;
    matrix_A=[matrix_A;full(sparse(1:sum(set_equal),find(set_equal),1))];
    [RREF,pivot_variables] = rref(matrix_A);
    free_variables =setdiff(1:size(matrix_A,2),pivot_variables);
    adj=sparse(size(RREF,2),size(free_variables,2));
    adj(pivot_variables,:)=-RREF(pivot_variables,free_variables);
    adj(free_variables,:)=speye(size(free_variables,2));
    mat_invhalf=sparse(size(free_variables,2),size(free_variables,2));
    set_nonzero_free=full(any(adj(pivot_variables,:),1));
    val=[];
    % nor_dis=[];
    % nor_dis2=[];
    t=1;
    sigma=1;
    sigma0=-Inf;
    count=0;
    for jj=0:max_iter % Max iter
        [f,g_x,hess]=func_all(x); %%%
        val=[val,f];      
        ve1=min(x-lb,ub-x);
        ve=1./(ve1.^2);
        mat=adj.'*spdiags(ve,0,size(ve1,1),size(ve1,1))*adj;   %%%
        try  %%%
            mat_invhalf(set_nonzero_free,set_nonzero_free)=inv(chol(adj(pivot_variables,set_nonzero_free).'*diag(ve(set_nonzero_free))*adj(pivot_variables,set_nonzero_free)+diag(ve(free_variables(set_nonzero_free)))));
            mat_invhalf(~set_nonzero_free,~set_nonzero_free)=spdiags(ve1(free_variables(~set_nonzero_free)),0,sum(~set_nonzero_free),sum(~set_nonzero_free));  %%%
        catch
            % keyboard
        end
        B=adj.'*hess*adj;
        g=adj.'*g_x;  %%%
        if norm(g)^2<opt_tol^2 % OptimalityTolerance
            fprintf('First optimality condition belows the tolerance. \n')
            break
        end
        s=mat_invhalf*cub_new_sub(mat_invhalf.'*g,mat_invhalf.'*B*mat_invhalf,sigma); %%%
        if ~isreal(s)
            % keyboard
        end
        if s.'*mat*s<thro
            delta_x=adj*s;
            count=count+1;
            sol_x=s;
        else
            sol_x_1=mynonconvexqp(adj.'*hess*adj,adj.'*g_x,mat,thro,mat_invhalf);   %%%
            try
            sol_x=sol_x_1;
            delta_x=adj*sol_x;
            catch
                keyboard
            end
        end
        if sol_x.'*mat*sol_x>thro+1e-1
            % keyboard
        end
        x1=x+t*delta_x;
        ff=func_obj(x); 
        f2=func_obj(x1);
        if any(x1<lb) || any (x1>ub)
            keyboard
            error('iteration point is not feasible \n')
        end
            
            
        %%%%%%%%%%%
        for jjj=1:100
            if f2>ff+1e-7%%%%+sol_x.'*g+0.5*sol_x.'*B*sol_x + sigma*(sol_x.'*mat*sol_x)^1.5/3 %%%%%%%
                sigma=2*sigma;
                s=mat_invhalf*cub_new_sub(mat_invhalf.'*g,mat_invhalf.'*B*mat_invhalf,sigma); %%%
                if s.'*mat*s<thro
                    delta_x=adj*s;
                    count=count+1;
                    sol_x=s;
                else
                    sol_x=sol_x_1;
                    delta_x=adj*sol_x;
                end
                x1=x+t*delta_x;
                f2=func_obj(x1);
            else
                break
            end
            if jjj==100 && f2>ff+1e-7
                % keyboard
            end
        end
        if ~isreal(delta_x)
            fprintf('The delta_x is not real. \n')
            keyboard
        end
        sigma=max(sigma/2,sigma0);
        x=x+t*delta_x;
        if norm(delta_x)<step_tol
            fprintf('Step size belows the tolerance. \n')
            break
        end
    end
end
