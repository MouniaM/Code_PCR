clc
clear all
%%
clc
clear all
% close all
DemoQatar='QatarDemo.xlsx';
PopSize= xlsread(DemoQatar);

InciQatar='DataCOVIDQatarHOA.xlsx';
IncidenceQatar= xlsread(InciQatar);

TestingDataQatarD='TestingDataQatar.xlsx';
TestingQatarData= xlsread(TestingDataQatarD);

%% Population size

nst=11; %% Number of status and stages
NAge=9; %% Number of Age groups
%% Initial population size and distribution
N0=zeros(NAge,1);
for a=2:2:16
    N0(round(a/2))=PopSize(a,1)+PopSize(a-1,1);  
end
N0(9)=PopSize(17,1);

%%%%%Reshaping the Data
mu=1/(80.7*365);               % 1/life expectancy in Qatar
eta=(ones(NAge,1)).*(1/(10*365)); %% Aging rate
eta(NAge)=0;
delta=1/(3.69);  % infection rate (Progression rate from exposed to infected)
nuM=1/(3.48); % Recovery Rate from mild infections
nuS=1/(4*7); %% Recovery rate from Severe disease
nuC=1/(6*7); %% Recovery rate from critical disease
nuSID=1/(3.48); %% Rate of progression from severe infection to severe disease
nuCID=1/(3.48); %% Rate of progression from critical infection to critical disease
DPCR=2; %% Prolonged PCR duration
nuP=1/(DPCR*7); 
DPAB=0.005; %% Duration of delay between recovery and antibody detection
nuPAB=1/(DPAB*7); 

load('EstimatedParamQatarNew.mat','EstimatedParamQatarNew') 

zz=EstimatedParamQatarNew;
epsiL=zz(2);
a0L=zz(6); 
Anch=zz(8);

FactorDeath=0.0002; %% Rate of mortality in those aged 30-39
Asymptomatic=0; %% it is the proportion of those who become asymptomatic; here is set to 0
Infectivity=1; %% Infectivity 
sigmaL=zeros(NAge,1); %% Susceptibility to infection from model fitiing qatar data
for j=1:NAge
    sigmaL(j)=zz(j+9);
end
fS1=0.011; %% Proportion of infections in those aged 30-39 that will progress to be infcetions that require hospitalizaton in acute care unit beds 
fC1=0.002; %% Proportion of infections in those aged 30-39 that will progress to be infcetions that require hospitalizaton in intensive care unit beds 
RRS=[0.1 0.1 0.5 1 1.2 2.3 4.5 7.8 27.6]; %%  Relative risk of severe infection
RRC=[0.2 0.2 0.3 1 1.8 4.7 10.6 13.6 8.7]; %%  Relative risk of critical infection
fS=zeros(NAge,1);
fC=zeros(NAge,1);
fM=zeros(NAge,1);
for a=1:NAge
fS(a)=RRS(a)*fS1; %% Proportion of infections in each age group that will progress to be infections that require hospitalizaton in acute care unit beds 
fC(a)=RRC(a)*fC1; %% Proportion of infections in each age group that will progress to be infections that require hospitalizaton in intensive care unit beds 
fM(a)=1-fS(a)-fC(a); %% Proportion of infections in each age group that will progress to be mild infections 
end

RRMor=[0.1 0.1 0.4 1 3 10 45 120 505]; %Relative Risk of morality
alpha=zeros(NAge,1); 
for a=1:NAge
    alpha(a)=FactorDeath*RRMor(a); %disease Mortality rate in each age group
end

%% Factor multiplying the contact rate and guiding to desired value of R0
ML=1.881 %% R0=3
% ML=0.99 %% R0=1.6;  
%%
%Initial conditions
initprev=1/9.*ones(NAge,1); %initial number of people with covid-19
x0L=zeros(NAge,nst);
for a=1:NAge
    x0L(a,1)= N0(a)-initprev(a);
    x0L(a,2)= 0+initprev(a);
end
    
x0=zeros(nst*NAge,1);
ii=1:nst;
for j=1:NAge
    isj= ii+(j-1)*nst;
    x0(isj)=x0L(j,ii);
end
%% Time scale (days)
t0=0;         % Start time
tf=2*365;         % Stop time (Based on Data Points)
dt=0.5;          % Time interval
tspan=t0:dt:tf;  % Timespan to use in the ode45
%% Resolution of the set of ODEs
[T,x] = ode45(@(t,x)risk_Corona19_PCR(t,x,dt,tspan,a0L,ML,sigmaL,alpha,delta,NAge,mu,epsiL,nst,eta,fM,fS,fC,nuM,nuS,nuC,nuSID,nuCID,nuP,nuPAB,Infectivity),tspan,x0);
TT=length(T);
for i=1:TT 
     [tx,lambdaLT]=risk_Corona19_PCR(T(i),x(i,:),dt,tspan,a0L,ML,sigmaL,alpha,delta,NAge,mu,epsiL,nst,eta,fM,fS,fC,nuM,nuS,nuC,nuSID,nuCID,nuP,nuPAB,Infectivity);
    for a=1:NAge
        lmbL(i,a)=lambdaLT(a);
    end
end
%% compartments by Age
SusA=x(:,1:nst:nst*NAge);                       %%Susceptible

LatentA=x(:,2:nst:nst*NAge);               %% Latently infected

InfectedIMA=x(:,3:nst:nst*NAge);  %% Mild infections

InfectedISA=x(:,4:nst:nst*NAge);  %% Severe Infections
InfectedDSA=x(:,5:nst:nst*NAge);  %% Severe Disease

InfectedICA=x(:,6:nst:nst*NAge);  %% Critical infections
InfectedDCA=x(:,7:nst:nst*NAge);  %% Critical Diseases

RecoveredA=x(:,8:nst:nst*NAge);   %% Recovered

PCRPA=x(:,9:nst:nst*NAge);        %% Those testing PCR positive after infectiousness
PABA=x(:,10:nst:nst*NAge);        %% Pre-antibody detection
ABA=x(:,11:nst:nst*NAge);         %% Antibody positive

TotalA=SusA+LatentA+InfectedIMA+InfectedISA+InfectedDSA++InfectedICA+InfectedDCA+RecoveredA;  %% Total Population
PropA=TotalA./(sum(TotalA,2)); %% proportion of each age group

InfA=(LatentA+InfectedIMA+InfectedISA+InfectedICA);               %% Infected in each age group
InfPCRA=(PCRPA+LatentA+InfectedIMA+InfectedISA+InfectedICA);      %% Infected and testing PCR positive after infectiousness in each age group

PropRecovA=RecoveredA./(TotalA);                %% Proportion of recovered in each age group
PropRecovABA=(ABA)./(TotalA);                   %% Porportion of testing antibody positive in each age group

PropRecov=sum(RecoveredA,2)./sum(TotalA,2);     %% Proportion of recovered in the total population
PropRecovAB=(sum(ABA,2))./sum(TotalA,2);        %% Porportion of testing antibody positive in the total population

Inf=sum(InfA,2);                                %% Infected in the total population
InfPCR=sum(InfPCRA,2);                          %% Infected and testing positive in the total population
PrevIA=(InfA)./TotalA;                          %% Prevalence of Infected in each age group                 
PrevIPCRA=(InfPCRA)./TotalA;                    %% Infected and testing PCR positive in each age group


PrevI=sum(InfA,2)./sum(TotalA,2);                       %% Prevalence of Infected in the total population 
PrevIPCR=sum(InfPCRA,2)./sum(TotalA,2);                 %% Prevalence of Infected and testing PCR positive after infectiousness in the total population 
for t=1:length(tspan)
    for a=1:NAge
   IncidenceAge(t,a)=delta*(LatentA(t,a));             %%Incidence in each age group
   CumIncidenceAge(t,a)=dt*sum(IncidenceAge(1:t,a));    %%CumulativeIncidence in each age group
   AttackRateTA(t,a)=100*CumIncidenceAge(t,a)/TotalA(t,a); %%Attack rate in each age group
    end
   Incidence(t)=sum(IncidenceAge(t,:));    %%Incidence in the total population
   AttackRateT(t)=100*sum(CumIncidenceAge(t,:))/sum(TotalA(t,:)); %%attack rate in the total population

end
%% Number of diagnosed cases :  true number of infected versus of actually observed number 
NbofTests=1000;
for t=1:length(tspan)
    NbIDiag(t)=binornd(NbofTests,PrevI(t)); %% True Number of diagnosed cases 
    NbPCRDiag(t)=binornd(NbofTests,PrevIPCR(t)); %% Actually-observed number of diagnosed cases
end
RPrev=PrevIPCR./PrevI;  % Ratio of prevalence 
RNbDiag=NbPCRDiag./NbIDiag; %% Ratio of diagnosed numbers


%% Number of diagnised number : True number of ever infected versus actually -observed number of seroprevalent
NbofTests=1000;
for t=1:length(tspan)
    NbRDiag(t)=binornd(NbofTests,PropRecov(t)); %% True number of ever infected
    NbABDiag(t)=binornd(NbofTests,PropRecovAB(t)); %% actually -observed number of seroprevalent
end
RPrev=PropRecovAB./PropRecov; %% Ratio of prevalence
RNbDiag=NbABDiag./NbRDiag; %% Ration of diagnosed cases


%% Figure Recovered/AB

%% Reproductive number
for a = 1:NAge
    omega(a)=0; %% omega in a rate at which an intervention is introduced here is set to 0
    beta(a)=ML.*a0L.*sigmaL(a);
   
   Term1_A(a,1) = (delta/(delta+mu+omega(a)+eta(a))).*(Asymptomatic*beta(a)/(nuM+mu+omega(a)+eta(a)));
   Term1_M(a,1) = (delta/(delta+mu+omega(a)+eta(a)))*((1-Asymptomatic)*fM(a)*beta(a)/(nuM+mu+omega(a)+eta(a)));
   Term1_S(a,1) = (delta/(delta+mu+omega(a)+eta(a)))*((1-Asymptomatic)*fS(a)*beta(a)/(nuSID+mu+omega(a)+eta(a)));
   Term1_C(a,1) = (delta/(delta+mu+omega(a)+eta(a)))*((1-Asymptomatic)*fC(a)*beta(a)/(nuCID+mu+omega(a)+eta(a)));
   
   Term2_A(a,1) = (omega(a)/(delta+mu+omega(a)+eta(a)))*(delta/(delta+mu+eta(a)))*(Asymptomatic*beta(a)/(nuM+mu+eta(a)));
   Term2_M(a,1) = (omega(a)/(delta+mu+omega(a)+eta(a)))*(delta/(delta+mu+eta(a)))*((1-Asymptomatic)*fM(a)*beta(a)/(nuM+mu+eta(a)));
   Term2_S(a,1) = (omega(a)/(delta+mu+omega(a)+eta(a)))*(delta/(delta+mu+eta(a)))*((1-Asymptomatic)*fS(a)*beta(a)/(nuSID+mu+eta(a)));
   Term2_C(a,1) = (omega(a)/(delta+mu+omega(a)+eta(a)))*(delta/(delta+mu+eta(a)))*((1-Asymptomatic)*fC(a)*beta(a)/(nuCID+mu+eta(a)));
    
   Term3_A(a,1) = (delta/(delta+mu+omega(a)+eta(a)))*(omega(a)/(nuM+mu+omega(a)+eta(a)))*(Asymptomatic*beta(a)/(nuM+mu+eta(a)));
   Term3_M(a,1) = (delta/(delta+mu+omega(a)+eta(a)))*(omega(a)/(nuM+mu+omega(a)+eta(a)))*((1-Asymptomatic)*fM(a)*beta(a)/(nuM+mu+eta(a)));
   Term3_S(a,1) = (delta/(delta+mu+omega(a)+eta(a)))*(omega(a)/(nuSID+mu+omega(a)+eta(a)))*((1-Asymptomatic)*fS(a)*beta(a)/(nuSID+mu+eta(a)));
   Term3_C(a,1) = (delta/(delta+mu+omega(a)+eta(a)))*(omega(a)/(nuCID+mu+omega(a)+eta(a)))*((1-Asymptomatic)*fC(a)*beta(a)/(nuCID+mu+eta(a)));
   
   R0A(a,1) = Term1_A(a,1)+Term1_M(a,1)+Term1_S(a,1)+Term1_C(a,1)+Term2_A(a,1)+Term2_M(a,1)+Term2_S(a,1)+Term2_C(a,1)+Term3_A(a,1)+Term3_M(a,1)+Term3_S(a,1)+Term3_C(a,1);
   R0_contribution(a,1) = PropA(1,a).*R0A(a,1);

end

R0= sum(R0_contribution)
