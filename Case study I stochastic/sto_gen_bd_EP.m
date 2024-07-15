%% Stochastic model generator
%  This is a data generator for stochastic model. 
%  Input:
%    - NR: (scalar) number of replicant
%    - D : (1 x ND vector) dosage levels
%    - T : (1 x NT vector) time points
%    - initial: (scalar) initial cell number
%    - Params: (1 x 6*S+1 vector) all necessary parameters ([p_i,r_i, d_i, b_i, E_i, n_i], c)
%    - num_sub: (scalar) number of subgroup
%
%  Output:
%

function DATA = sto_gen_bd_EP(NR,Conc,Time,initial,theta)
    %% Initialize the parameter
    c      = theta(end);
    theta(end) = [];
    PMAT   = reshape(theta,6,[])'; % S x 6 parameters matrix
    p_vec  = PMAT(:,1);
    PMAT(:,1) = [];
    ND     = length(Conc);
    NT     = length(Time);
    DATA    = zeros(NT,ND,NR);
    for i = 1:NR
        for j = 1:ND
            DATA(1,j,i) = initial;
            for k = 2:NT
                temp_time = [Time(1),Time(k)];
                temp = Sto_samplepath(initial,PMAT',Conc(j),temp_time,p_vec');
                DATA(k,j,i) = temp(2);
            end
        end
    end
    DATA(2:end,:,:) = max(DATA(2:end,:,:) + round(normrnd(0,c,NT-1,ND,NR)),0);
end