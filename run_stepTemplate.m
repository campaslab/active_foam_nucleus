% function ns_simDynamics(bs,delta,psi,mu,taut,nvf,...
%     nstf,taun,taua,mun,nrst,asp,runTm,tmSt,svPt,edpc,rpc)


%% Initialize workspace
clearvars;      close all;

%% Run simulations

bs=6;
delta=10;
psi=1;
mu=0.75;
taut=10;

taun=1;
taua=1;
mun=0;
nrst=1;
nstf=50;
asp=1;
runTm=60;
tmSt=0.01;
svPt=10;
edpc=0.15;
ststn=1.5;

for rpc=2:5
    tic
    nvf=0.85;
    
    ns_simStep(bs,delta,psi,mu,taut,nvf,...
        nstf,taun,taua,mun,nrst,asp,ststn,runTm,tmSt,svPt,edpc,rpc);
    toc
    
%     nvf=0;
% 
%     tic
%     ns_simStep(bs,delta,psi,mu,taut,nvf,...
%         nstf,taun,taua,mun,nrst,asp,ststn,runTm,tmSt,svPt,edpc,rpc);
%     toc
end