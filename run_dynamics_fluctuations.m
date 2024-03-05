% function ns_simDynamics(bs,delta,psi,mu,taut,nvf,...
%     nstf,taun,taua,mun,nrst,asp,runTm,tmSt,svPt,edpc,rpc)


%% Initialize workspace
clearvars;      close all;

tautv=[1,2,5,10];

for tac=2:4
%% Run simulations

    [bs,delta,psi,mu,taut,nvf,nstf,taun,taua,mun,nrst,asp,runTm,tmSt,svPt,edpc,rpc]=...
        deal(6,10,1,0.05*tautv(tac),tautv(tac),0,0,1,1,0,1,1,30,0.01,10,0.15,1);

    tic
    ns_simDynamics(bs,delta,psi,mu,taut,nvf,...
        nstf,taun,taua,mun,nrst,asp,runTm,tmSt,svPt,edpc,rpc);
    toc
end