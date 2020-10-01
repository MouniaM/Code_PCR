function [Lambda]=ForceofInfection(x,betaL,epsiL,Nage,Infectivity)

eps=1e-10;

%Force of infection
%rhoN=zeros(Nage,1);
rhoN=sum(x(:,1:8),2); %%Population size in each age group
rhoNtot=sum(rhoN(1:Nage));

%Mixing Matrix
%%%Age group
II=eye(Nage);
H=zeros(Nage,Nage);
for a=1:Nage
    H(:,a)=epsiL.*II(:,a)+(1-epsiL).*(rhoN(a)/(rhoNtot+eps)).*ones(Nage,1); %% Age mixing matrix
end

LambdaR=zeros(Nage,1);

for a=1:Nage %Age
    LambdaR(:)=LambdaR(:)+betaL(:).*H(:,a).*((x(a,3)+x(a,4)+x(a,6))/(rhoN(a)+eps)); %% Force of infection
end

Lambda=LambdaR; 


