function NEWTREAT = generateData_param(param,Conc,Time,NR,num_sub,init)

[Pop_vec, Param_mat, noise] = reshape_param(param,num_sub);


NC = length(Conc);
NT = length(Time);
NEWTREAT = zeros(NR, NC, NT); %store cell counts (NR is number of replicates)

MEAN_initial = init;
for j=1:NR %replicate
%      eps =0;
     for k=1:NT %time
        for n=1:NC %concentration
            
            for ii = 1:num_sub
                NEWTREAT(j,n,k) = NEWTREAT(j,n,k)+MEAN_initial*Pop_vec(ii)*popfunc(Param_mat(:,ii),Conc(n),Time(k));
%                 NEWTREAT(j,n,k) = (MEAN_initial + eps)*Pop1*popfunc(min_beta_E_827,Conc(n),Time(k)) + ...
%                     (MEAN_initial)*Pop2*popfunc(min_beta_E_1975,Conc(n),Time(k))+...
%                     (MEAN_initial)*Pop3*popfunc(min_beta_P_1975,Conc(n),Time(k))+...
%                     (MEAN_initial)*(1-Pop1-Pop2-Pop3)*popfunc(Pop4Params,Conc(n),Time(k));
            end
            if ((Conc(n)<=0.1) && (Time(k)>=24))
                
                %NEWTREAT(j,n,k) = NEWTREAT(j,n,k)*(1 + (rand*NoiseH - NoiseH/2));
                NEWTREAT(j,n,k) = max(NEWTREAT(j,n,k) + (randn*max(noise)), 0);
                
            else
                 %NEWTREAT(j,n,k) = NEWTREAT(j,n,k)*(1 + (rand*NoiseL - NoiseL/2));
                NEWTREAT(j,n,k) = max(NEWTREAT(j,n,k)+ (randn*min(noise)),0);
            end
        end

    end
end





end