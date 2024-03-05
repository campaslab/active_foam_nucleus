function ns_simStep(bs,delta,psi,mu,taut,nvf,...
    nstf,taun,taua,mun,nrst,asp,ststn,runTm,tmSt,svPt,edpc,rpc)

% clearvars;
% [bs,delta,psi,mu,taut,nvf,nstf,taun,taua,mun,nrst,asp,runTm,tmSt,svPt,edpc,rpc]=...
%     deal(6,10,1,1.25,10,0.8,50,1,1,0,1,0.25,40,0.01,10,0.15,1);

%% Specify parameters for initical configurations
% Here, we categorize these parameterse into two categories, geometric
% parameters and mechanical parameters. Denote a set of geometric
% parameters as 'gm_p' and a set of mechanical parameters as 'mc_p'. If
% we want to add paramters, just add more elements in each variable.

[gmp,mcp]=deal(struct);
[gmp.tgAr,gmp.sstn,gmp.bs,gmp.edpc,gmp.asp]=deal(1,0,bs,edpc,asp);
[mcp.del,mcp.psi,mcp.mu,mcp.taut,mcp.nvf,mcp.nstf,mcp.taun,...
    mcp.taua,mcp.mun,mcp.nrst]=...
    deal(delta,psi,mu,taut,nvf,nstf,taun,taua,mun,nrst);
[gmp.nFa,gmp.lsc,gmp.dt]=...
    deal(gmp.bs^2,mcp.psi*sqrt(gmp.tgAr),tmSt);
[gmp.major,gmp.minor]=deal(sqrt(mcp.nvf*mcp.psi/pi/gmp.asp),...
    sqrt(mcp.nvf*mcp.psi/pi*gmp.asp));
[gmp.runTm,gmp.svPt]=deal(runTm,svPt);

gmp.shEd=0.01*2*sqrt(pi);

gmp.ellDis=[gmp.major*(cos(0:pi/30:2*pi-pi/30)).',...
    gmp.minor*(sin(0:pi/30:2*pi-pi/30)).';gmp.major,0];
gmp.ellAng=[mod(atan2(gmp.ellDis(1:end-1,2),gmp.ellDis(1:end-1,1)),2*pi);2*pi];
gmp.ellDis=sqrt(sum(gmp.ellDis.^2,2));

%% Generate initial configuration
[vrtx,edge,face,rg]=ns_initialConf(gmp.bs); 

%% Relaxation of the configuration with extracellular spaces. 
[itc,ttMx,ttTm]=deal(0,10*mcp.taut,0);

gmp.eLnc=ns_edgeLenAll(edge,gmp);

while (ttTm<=ttMx)    
    [vrtx,edge,face]=ns_iteration(vrtx,edge,face,rg,gmp,mcp);
    itc=itc+1;
    ttTm=ttTm+gmp.dt;
    
    if mod(itc,5)==0
        eLn=ns_edgeLenAll(edge,gmp);
        t1Id=find(eLn<gmp.shEd);
        if isempty(t1Id)==0
            for jj=1:size(t1Id,1)
                chk=min(face{1}(edge{1}(t1Id(jj),rg.ei(2):rg.ef(2)),rg.fi(5)));
                if chk>3 && gmp.eLnc(t1Id(jj))>eLn(t1Id(jj))
                    [vrtx,edge,face]=...
                        ns_t1Transition(vrtx,edge,face,rg,gmp,mcp,t1Id(jj));
                end
            end
        end
        gmp.eLnc=ns_edgeLenAll(edge,gmp);
    end       
    
    %% Compute t4t1 transition
    for ii=1:size(vrtx,1)
        [vrtx,edge,face]=...
            ns_t4AdjacentTransition(vrtx,edge,face,rg,gmp,mcp,ii);
    end
end
% 
%% Relaxation of the configuration with extracellular spaces. 
[itc,ttMx,ttTm,exc]=deal(0,runTm*mcp.taut,0,1);    
ttDtPt=ttMx/gmp.dt/svPt;
dtPt=4000;
dtc=1;    

[vr,ed1,ed2,ed3,ed4,fa1,fa2,fa3]=deal(cell(min(dtPt,ttDtPt),1));
stnAll=zeros(min(dtPt,ttDtPt),1);

odir=sprintf('step/nvf%s_nstf%s_mu%s_mun%s_taun%s_taut%s_nrst%s_asp%s/',num2str(nvf),...
    num2str(nstf),num2str(mu),num2str(mun),num2str(taun),num2str(taut),num2str(nrst),num2str(asp));
odir=strrep(odir,'.','d');
if ~exist(odir,'dir')
    mkdir(odir);
end


while (ttTm<=2*mcp.taut)
    [vrtx,edge,face]=ns_iteration(vrtx,edge,face,rg,gmp,mcp);
    itc=itc+1;
    ttTm=ttTm+gmp.dt;
    
    if mod(itc,5)==0
        eLn=ns_edgeLenAll(edge,gmp);
        t1Id=find(eLn<gmp.shEd);
        if isempty(t1Id)==0
            for jj=1:size(t1Id,1)
                chk=min(face{1}(edge{1}(t1Id(jj),rg.ei(2):rg.ef(2)),rg.fi(5)));
                if chk>3 && gmp.eLnc(t1Id(jj))>eLn(t1Id(jj))
                    [vrtx,edge,face]=...
                        ns_t1Transition(vrtx,edge,face,rg,gmp,mcp,t1Id(jj));
                end
            end
        end
        gmp.eLnc=ns_edgeLenAll(edge,gmp);
    end
    
    if mod(itc,5)==2
        %% Compute t4t1 transition
        for ii=1:size(vrtx,1)
            [vrtx,edge,face]=...
                ns_t4AdjacentTransition(vrtx,edge,face,rg,gmp,mcp,ii);
        end
    end

    if mod(itc,svPt)==3
        
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

face{1}(:,rg.fi(7))=face{1}(:,rg.fi(7))+ststn*face{1}(:,rg.fi(7)+1);

%% Relax configuration
while (ttTm<=ttMx)
    [vrtx,edge,face]=ns_iteration(vrtx,edge,face,rg,gmp,mcp);
    itc=itc+1;
    ttTm=ttTm+gmp.dt;
    
    if mod(itc,5)==0
        eLn=ns_edgeLenAll(edge,gmp);
        t1Id=find(eLn<gmp.shEd);
        if isempty(t1Id)==0
            for jj=1:size(t1Id,1)
                chk=min(face{1}(edge{1}(t1Id(jj),rg.ei(2):rg.ef(2)),rg.fi(5)));
                if chk>3 && gmp.eLnc(t1Id(jj))>eLn(t1Id(jj))
                    [vrtx,edge,face]=...
                        ns_t1Transition(vrtx,edge,face,rg,gmp,mcp,t1Id(jj));
                end
            end
        end
        gmp.eLnc=ns_edgeLenAll(edge,gmp);
    end
    
    if mod(itc,5)==2
        %% Compute t4t1 transition
        for ii=1:size(vrtx,1)
            [vrtx,edge,face]=...
                ns_t4AdjacentTransition(vrtx,edge,face,rg,gmp,mcp,ii);
        end
    end

    if mod(itc,svPt)==3
        
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
    vr=vr(1:exc-1);
    ed1=ed1(1:exc-1);
    ed2=ed2(1:exc-1);
    ed3=ed3(1:exc-1);
    ed4=ed4(1:exc-1);
    fa1=fa1(1:exc-1);
    fa2=fa2(1:exc-1);
    fa3=fa3(1:exc-1);
    
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

ns_exportParameter(gmp,mcp,odir,rpc);
% 
end
    