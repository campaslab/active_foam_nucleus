%% Step simulations
function ts_simStep(bs,beta,gamma,delta,psi,mu,taut,edpc,rpc,ststn)

%% Specify parameters for initical configurations
% Here, we categorize these parameterse into two categories, geometric
% parameters and mechanical parameters. Denote a set of geometric
% parameters as 'gm_p' and a set of mechanical parameters as 'mc_p'. If
% we want to add paramters, just add more elements in each variable.

[gmp,mcp]=deal(struct);
[gmp.tgAr,gmp.sstn,gmp.bs,gmp.edpc]=deal(1,0,bs,edpc);
[mcp.beta,mcp.gam,mcp.del,mcp.psi,mcp.mu,mcp.taut]=...
    deal(beta,gamma,delta,psi,mu,taut);
[gmp.nFa,gmp.lsc,gmp.dt]=...
    deal(gmp.bs^2,mcp.psi*sqrt(gmp.tgAr),0.005);
gmp.shEd=0.01*2*sqrt(pi);


%% Generate initial configuration
[vrtx,edge,face,rg]=ts_initialConf(gmp.bs);

% set an initial tension as fixed point tension
for ii=1:size(edge{1},1)
    edge{1}(ii,rg.ei(3))=ts_edgeFixedTension(face,edge,rg,mcp,ii);
end

%% Introduce extracellular space cells
for ii=1:size(vrtx,1)
    [vrtx,edge,face,gmp]=...
        ts_t2ReverseTransition(vrtx,edge,face,rg,gmp,mcp,ii);              
end

%% Relaxation of the configuration with extracellular spaces. 
mcp.mu=mu;
[itc,ttMx,ttTm]=deal(0,4*mcp.taut+8*gmp.dt,0);
gmp.shEd=0.01*2*sqrt(pi);

while (ttTm<=ttMx)    
    gmp.eLnc=ts_edgeLenAll(edge,gmp);

    [vrtx,edge,face]=ts_iteration(vrtx,edge,face,rg,gmp,mcp);
    itc=itc+1;
    ttTm=ttTm+gmp.dt;

    if mod(itc,3)==0
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType1(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,5)==1 
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType3(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,5)==2
        %% Compute t4t1 transition
        for ii=1:size(vrtx,1)
            [vrtx,edge,face]=...
                ts_t4AdjacentTransition(vrtx,edge,face,rg,gmp,mcp,ii);
        end
    end

    if mod(itc,5)==3 
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType2(vrtx,edge,face,rg,gmp,mcp);
    end
end

%% Relaxation of the configuration with extracellular spaces. 
[itc,ttMx,ttTm,exc]=deal(0,60*mcp.taut,0,1);    
ttDtPt=ttMx/gmp.dt/10;
dtPt=4000;
dtc=1;    

[vr,ed1,ed2,ed3,ed4,fa1,fa2,fa3]=deal(cell(min(dtPt,ttDtPt),1));
stnAll=zeros(min(dtPt,ttDtPt),1);

odir=sprintf('step/psi%s_gam%s_mu%s/',num2str(psi),...
    num2str(gamma),num2str(mu));
odir=strrep(odir,'.','d');
if ~exist(odir,'dir')
    mkdir(odir);
end

while (ttTm<=2*mcp.taut)
    gmp.eLnc=ts_edgeLenAll(edge,gmp);

    [vrtx,edge,face]=ts_iteration(vrtx,edge,face,rg,gmp,mcp);
    itc=itc+1;
    ttTm=ttTm+gmp.dt;

    if mod(itc,3)==0
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType1(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,5)==1 
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType3(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,5)==2
        %% Compute t4t1 transition
        for ii=1:size(vrtx,1)
            [vrtx,edge,face]=...
                ts_t4AdjacentTransition(vrtx,edge,face,rg,gmp,mcp,ii);
        end
    end

    if mod(itc,5)==3 
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType2(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,10)==3
        [vr{exc},ed1{exc},ed2{exc},ed3{exc},ed4{exc},fa1{exc},...
            fa2{exc},fa3{exc}]=deal(vrtx,edge{1},edge{2},edge{3},...
            edge{4},face{1},face{2},face{3});
        stnAll(exc)=gmp.sstn;
        exc=exc+1;
        if exc==dtPt+1
            save(sprintf('%s/vr_r%d_d%d.mat',odir,rpc,dtc),'vr');
            save(sprintf('%s/ed1_r%d_d%d.mat',odir,rpc,dtc),'ed1');
            save(sprintf('%s/ed2_r%d_d%d.mat',odir,rpc,dtc),'ed2');
            save(sprintf('%s/ed3_r%d_d%d.mat',odir,rpc,dtc),'ed3');
            save(sprintf('%s/ed4_r%d_d%d.mat',odir,rpc,dtc),'ed4');
            save(sprintf('%s/fa1_r%d_d%d.mat',odir,rpc,dtc),'fa1');
            save(sprintf('%s/fa2_r%d_d%d.mat',odir,rpc,dtc),'fa2');
            save(sprintf('%s/fa3_r%d_d%d.mat',odir,rpc,dtc),'fa3');
            save(sprintf('%s/stn_r%d_d%d.mat',odir,rpc,dtc),'stnAll');
            dtc=dtc+1;
            exc=1;
            ttDtPt=ttDtPt-dtPt;
            [vr,ed1,ed2,ed3,ed4,fa1,fa2,fa3]=...
                deal(cell(min(dtPt,ttDtPt),1)); 
            stnAll=zeros(min(dtPt,ttDtPt),1);
        end    
    end
end

%% Apply step strain
gmp.sstn=ststn;
vrtx(:,rg.vi(3))=vrtx(:,rg.vi(3))+ststn*vrtx(:,rg.vi(3)+1);

for edc=1:size(edge{2},1)
    emd=edge{2}{edc};
    emd(:,1)=emd(:,1)+ststn*emd(:,2);
    edge{2}{edc}=emd;
end

%% Relax configurations.    
while (ttTm<=ttMx)
    gmp.eLnc=ts_edgeLenAll(edge,gmp);

    [vrtx,edge,face]=ts_iteration(vrtx,edge,face,rg,gmp,mcp);
    itc=itc+1;
    ttTm=ttTm+gmp.dt;

    if mod(itc,3)==0
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType1(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,5)==1 
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType1(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,5)==2
        %% Compute t4t1 transition
        for ii=1:size(vrtx,1)
            [vrtx,edge,face]=...
                ts_t4AdjacentTransition(vrtx,edge,face,rg,gmp,mcp,ii);
        end
    end

    if mod(itc,5)==3 
        [vrtx,edge,face,gmp]=...
            ts_topTransitionType1(vrtx,edge,face,rg,gmp,mcp);
    end

    if mod(itc,10)==3
        [vr{exc},ed1{exc},ed2{exc},ed3{exc},ed4{exc},fa1{exc},...
            fa2{exc},fa3{exc}]=deal(vrtx,edge{1},edge{2},edge{3},...
            edge{4},face{1},face{2},face{3});
        stnAll(exc)=gmp.sstn;
        exc=exc+1;
        if exc==dtPt+1
            save(sprintf('%s/vr_r%d_d%d.mat',odir,rpc,dtc),'vr');
            save(sprintf('%s/ed1_r%d_d%d.mat',odir,rpc,dtc),'ed1');
            save(sprintf('%s/ed2_r%d_d%d.mat',odir,rpc,dtc),'ed2');
            save(sprintf('%s/ed3_r%d_d%d.mat',odir,rpc,dtc),'ed3');
            save(sprintf('%s/ed4_r%d_d%d.mat',odir,rpc,dtc),'ed4');
            save(sprintf('%s/fa1_r%d_d%d.mat',odir,rpc,dtc),'fa1');
            save(sprintf('%s/fa2_r%d_d%d.mat',odir,rpc,dtc),'fa2');
            save(sprintf('%s/fa3_r%d_d%d.mat',odir,rpc,dtc),'fa3');
            save(sprintf('%s/stn_r%d_d%d.mat',odir,rpc,dtc),'stnAll');
            dtc=dtc+1;
            exc=1;
            ttDtPt=ttDtPt-dtPt;
            [vr,ed1,ed2,ed3,ed4,fa1,fa2,fa3]=...
                deal(cell(min(dtPt,ttDtPt),1)); 
            stnAll=zeros(min(dtPt,ttDtPt),1);
        end    
    end
end

if exc>1
    save(sprintf('%s/vr_r%d_d%d.mat',odir,rpc,dtc),'vr');
    save(sprintf('%s/ed1_r%d_d%d.mat',odir,rpc,dtc),'ed1');
    save(sprintf('%s/ed2_r%d_d%d.mat',odir,rpc,dtc),'ed2');
    save(sprintf('%s/ed3_r%d_d%d.mat',odir,rpc,dtc),'ed3');
    save(sprintf('%s/ed4_r%d_d%d.mat',odir,rpc,dtc),'ed4');
    save(sprintf('%s/fa1_r%d_d%d.mat',odir,rpc,dtc),'fa1');
    save(sprintf('%s/fa2_r%d_d%d.mat',odir,rpc,dtc),'fa2');
    save(sprintf('%s/fa3_r%d_d%d.mat',odir,rpc,dtc),'fa3');
    save(sprintf('%s/stn_r%d_d%d.mat',odir,rpc,dtc),'stnAll');
end    

end