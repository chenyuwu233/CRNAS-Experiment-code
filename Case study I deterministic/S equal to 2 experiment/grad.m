%% Get the total grad

function [g_p,g_param,hess]=grad(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub,opt)
% val=MEAN_initial*Pop(1)*popfunc(param[:,1],Conc,Time.')+ ...
%                 MEAN_initial*Pop(2)*popfunc(param[:,2],Conc,Time')+...
%                 MEAN_initial*Pop(3)*popfunc(param[:,3],Conc,Time')+...
%                 MEAN_initial*Pop(4)*popfunc(param[:,4],Conc,Time');
% val_1=MEAN_initial*popfunc(param[:,1],Conc,Time.');
% g_1=mean(mean((val-Observations_ave).*val_1));
% g=(val-Observations_ave)*MEAN_initial*Pop(1).*val_1.*Time.';
% de=(Conc./param[3,1]).^param[4,1];
% g=(val-Observations_ave)*MEAN_initial*Pop(1).*val_1.*Time.'.*de./(param[2,1].*de+1);

if nargin <9
    for ii=1:num_sub
        [pop(:,:,ii),g(:,:,1,ii),g(:,:,2,ii),g(:,:,3,ii),g(:,:,4,ii)] = grad_pop(param(:,ii),Conc,Time.');
    end
    %val=MEAN_initial*sum(permute(permute(pop,[3,1,2]).*Pop, [2,3,1]),3);
    Pop3(1,1,:)=Pop;
    val=MEAN_initial*sum(pop.*Pop3,3);
    g_p=squeeze(mean(MEAN_initial*pop.*(val-Observations_ave),[1,2]));
    g_param=squeeze(mean(MEAN_initial*g.*(val-Observations_ave),[1,2]));
    g_param=g_param.*Pop.';
else
    for ii=1:num_sub
        [pop(:,:,ii),g(:,:,1,ii),g(:,:,2,ii),g(:,:,3,ii),g(:,:,4,ii),hess(:,:,:,:,ii)] = grad_pop(param(:,ii),Conc,Time.',2);
    end
    Pop3(1,1,:)=Pop;
    val=MEAN_initial*sum(pop.*Pop3,3);
    g_p=squeeze(mean(MEAN_initial*pop.*(val-Observations_ave),[1,2]));
    g_param=squeeze(mean(MEAN_initial*g.*(val-Observations_ave),[1,2]));
    g_param=g_param.*Pop.';
    
    
    hess_pp(:,:,1,:)=pop;
    hess_pp=hess_pp.*pop;
    hess_pp=squeeze(MEAN_initial^2*mean(hess_pp,[1,2]));
    
    err=Observations_ave-val;
    Pop4(1,1,1,:)=Pop;
    gp2=permute(reshape(permute(g.*Pop4,[3,4,1,2]),[],size(g,1),size(g,2)),[2,3,1]);
    gp3(:,:,1,:)=gp2;
    hess_ppar=squeeze(MEAN_initial^2*mean(pop.*gp3,[1,2]));
    if size(hess_ppar,1)~= num_sub
        hess_ppar = hess_ppar';
    end

    for ii=1:num_sub
        hess_ppar(ii,(4*ii-3):4*ii) = hess_ppar(ii,(4*ii-3):4*ii)-MEAN_initial*squeeze(mean(err.*g(:,:,:,ii),[1,2])).';
    end
    


    hess_parpar = squeeze(MEAN_initial^2*mean(gp2.*gp3,[1,2]));
    for ii=1:num_sub
        hess_parpar((4*ii-3):4*ii,(4*ii-3):4*ii) = hess_parpar((4*ii-3):4*ii,(4*ii-3):4*ii)-MEAN_initial*Pop(ii)*squeeze(mean(err.*hess(:,:,:,:,ii),[1,2]));
    end

    hess=[hess_pp,hess_ppar;hess_ppar.',hess_parpar];
    
end

end