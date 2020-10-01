function [dxL,lambda,lambdaNLT]=risk_Corona19_PCR(t,x,dt,tspan,a0L,ML,sigmaL,alpha,delta,NAge,mu,epsiL,nst,eta,fM,fS,fC,nuM,nuS,nuC,nuSID,nuCID,nuP,nuPAB,Infectivity)
%% %%%%%%%%Gradually Testing
xt=zeros(NAge,nst);
ii=1:nst;
for j=1:NAge
    isj= ii+(j-1)*nst;
    xt(j,ii)=x(isj);
end


%% %% %%%%%%%%%% Contact rate
betaL=zeros(NAge,1);
for a=1:NAge
 betaL(a)=ML.*a0L.*sigmaL(a);
end


%% 
lambda=ForceofInfection(xt,betaL,epsiL,NAge,Infectivity);
%%
%
%Age 0-4
dx=zeros(NAge,nst);
%%%%%%%%People in Lock-Down
dx(1,1)=-(mu+lambda(1)+eta(1))*xt(1,1);                                                              %(S)

dx(1,2)= lambda(1)*xt(1,1) - (mu+delta+eta(1))*xt(1,2);                                              %(E)

dx(1,3)=fM(1)*delta*xt(1,2) - (mu+nuM+eta(1))*xt(1,3);                                               %(IM)

dx(1,4)=fS(1)*delta*xt(1,2) - (mu+nuSID+eta(1))*xt(1,4);                                             %(IS) 
dx(1,5)=nuSID*xt(1,4) - (mu+nuS+eta(1))*xt(1,5);                                                     %(DS)

dx(1,6)=fC(1)*delta*xt(1,2) - (mu+nuCID+eta(1))*xt(1,6);                                             %(IC)
dx(1,7)=nuCID*xt(1,6) - (mu+nuC+eta(1)+alpha(1))*xt(1,7);                                            %(DC)

dx(1,8)=nuM*xt(1,3)+nuS*xt(1,5)+nuC*xt(1,7)-(mu+eta(1))*xt(1,8);                                     %(RSym)
dx(1,9)=nuM*xt(1,3)+nuSID*xt(1,4)+nuCID*xt(1,6)-(mu+nuP+eta(1))*xt(1,9)-alpha(1)*xt(1,7);            %(YPCR)
dx(1,10)=nuM*xt(1,3)+nuSID*xt(1,4)+nuCID*xt(1,6)-(mu+nuPAB+eta(1))*xt(1,10)-alpha(1)*xt(1,7);        %(XIgG)
dx(1,11)=nuPAB*xt(1,10)-(mu+eta(1)).*xt(1,11);                                                       %(ZIgG)


for rst=2:NAge
    %%%%%%%%People in Lock-Down
dx(rst,1)=eta(rst-1).*xt(rst,1)- (mu+lambda(rst)+eta(rst))*xt(rst,1);                                                                          %(S)

dx(rst,2)=eta(rst-1).*xt(rst,2)+lambda(rst)*xt(rst,1) - (mu+delta+eta(rst))*xt(rst,2) ;                                                        %(E)

dx(rst,3)=eta(rst-1).*xt(rst,3)+fM(rst)*delta*xt(rst,2) - (mu+nuM+eta(rst))*xt(rst,3);                                                         %(IM)

dx(rst,4)=eta(rst-1).*xt(rst,4)+fS(rst)*delta*xt(rst,2) - (mu+nuSID+eta(rst))*xt(rst,4);                                                       %(IS) 
dx(rst,5)=eta(rst-1).*xt(rst,5)+nuSID*xt(rst,4) - (mu+nuS+eta(rst))*xt(rst,5);                                                                 %(DS)

dx(rst,6)=eta(rst-1).*xt(rst,6)+fC(rst)*delta*xt(rst,2) - (mu+nuCID+eta(rst))*xt(rst,6);                                                       %(IC)
dx(rst,7)=eta(rst-1).*xt(rst,7)+nuCID*xt(rst,6) - (mu+nuC+eta(rst)+alpha(rst))*xt(rst,7);                                                      %(DC)

dx(rst,8)=eta(rst-1).*xt(rst,8)+nuM*xt(rst,3)+nuS*xt(rst,5)+nuC*xt(rst,7)-(mu+eta(rst))*xt(rst,8);                                             %(RSym)
dx(rst,9)=eta(rst-1).*xt(rst,9)+nuM*xt(rst,3)+nuSID*xt(rst,4)+nuCID*xt(rst,6)-(mu+nuP+eta(rst))*xt(rst,9)-alpha(rst)*xt(rst,7);                %(YPCR)
dx(rst,10)=eta(rst-1).*xt(rst,10)+nuM*xt(rst,3)+nuSID*xt(rst,4)+nuCID*xt(rst,6)-(mu+nuPAB+eta(rst))*xt(rst,10)-alpha(rst)*xt(rst,7);           %(XIgG)
dx(rst,11)=eta(rst-1).*xt(rst,11)+nuPAB*xt(rst,10)-(mu+eta(rst)).*xt(rst,11);                                                                  %(ZIgG)
end


dxL=zeros(NAge*nst,1);
ii=1:nst;
for j=1:NAge
    isj=ii+(j-1)*nst;
dxL(isj)= dx(j,ii);
end

end



