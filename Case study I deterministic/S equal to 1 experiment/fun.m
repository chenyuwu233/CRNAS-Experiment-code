%% Obtain the total function value

function f=fun(Conc,Time,NR,Pop,param,Observations_ave,MEAN_initial,num_sub)
%val=MEAN_initial*Pop(1)*popfunc(param(:,1),Conc,Time.')+ ...
%                 MEAN_initial*Pop(2)*popfunc(param(:,2),Conc,Time')+...
%                 MEAN_initial*Pop(3)*popfunc(param(:,3),Conc,Time')+...
%                 MEAN_initial*Pop(4)*popfunc(param(:,4),Conc,Time');

for ii=1:num_sub
pop(:,:,ii) = popfunc2(param(:,ii),Conc,Time.');
end
val=MEAN_initial*sum(permute(permute(pop,[3,1,2]).*Pop, [2,3,1]),3);

f=mean((val-Observations_ave).^2,'all')/2;
end