%% Adjust vertex position in force direction.
% For a given set of vertex coordinates, calculate force for each vertex
% and displace vertex by force*scFac*edLnAv. 

function [vrtxN,edgeN,faceN]=ns_iteration(vrtx,edge,face,rg,gmp,mcp)

% tic;

%% Declare new variables
[vrtxN,edgeN,faceN]=deal(vrtx,edge,face);
[vMx,eMx]=deal(size(vrtx,1),size(edge{1},1));
vrtxN(:,rg.vi(4):rg.vf(4))=zeros(vMx,2);
faceN{1}(:,rg.fi(8):rg.ff(8))=zeros(gmp.nFa,2);
faceN{1}(:,rg.fi(10):rg.ff(10))=zeros(gmp.nFa,1);

%% Compute net force for intermediate vertices
for edc=1:eMx
    rfn=edge{1}(edc,rg.ei(4));
    vTn=zeros(rfn+2,2);
    eMd=edge{2}{edc};
    fCn=face{1}(edge{1}(edc,4),rg.fi(6):rg.ff(6));
    eTf=max(0,edge{1}(edc,rg.ei(3)));

    for jj=1:rfn+1
        [tv,nv,eln]=ns_edgeVector(eMd(jj,:),eMd(jj+1,:),...
            fCn,gmp.bs,gmp.sstn);
        eNf=ns_edgeNormalForce(edge,face,rg,edc,mcp.del,eln,gmp.lsc);
        vTn(jj,:)=vTn(jj,:)+eTf*tv;
        vTn(jj+1,:)=vTn(jj+1,:)-eTf*tv;
        vTn(jj,:)=vTn(jj,:)+eNf*nv;
        vTn(jj+1,:)=vTn(jj+1,:)+eNf*nv;
    end
    vrtxN(edgeN{1}(edc,1),rg.vi(4):rg.vf(4))=...
        vrtxN(edgeN{1}(edc,1),rg.vi(4):rg.vf(4))+vTn(1,:);
    vrtxN(edgeN{1}(edc,2),rg.vi(4):rg.vf(4))=...
        vrtxN(edgeN{1}(edc,2),rg.vi(4):rg.vf(4))+vTn(end,:);

    edgeN{3}{edc}=vTn;        
    edgeN{2}{edc}=eMd+mcp.psi*gmp.dt*vTn;
end

%% Compute nucleus forces for each edge

if mcp.nvf>0    
    for fac=1:gmp.nFa
        feId=abs(face{3}{fac});
        fcen=face{1}(fac,rg.fi(7):rg.ff(7));
        pcen=face{1}(fac,rg.fi(6):rg.ff(6));
        forn=face{1}(fac,rg.fi(9));
        for edc=1:size(feId,2)
            eMd=edge{2}{feId(edc)};
            ecDis=[fcen;pcen;eMd];
            ecDis=ns_crdLocal(ecDis,gmp.bs,gmp.sstn);        
            fcen=ecDis(1,:);
            pcen=ecDis(2,:);
            eMd=ecDis(3:end,:);
            ecDis=sqrt(sum((eMd-fcen).^2,2));

%             if min(ecDis)<eqDis
            [eFrc,nFrc,ncMmt]=ns_nucleusForce(gmp,mcp,fcen,pcen,forn,eMd,ecDis);
            vrtxN(edgeN{1}(feId(edc),1),rg.vi(4):rg.vf(4))=...
                vrtxN(edgeN{1}(feId(edc),1),rg.vi(4):rg.vf(4))+eFrc(1,:);
            vrtxN(edgeN{1}(feId(edc),2),rg.vi(4):rg.vf(4))=...
                vrtxN(edgeN{1}(feId(edc),2),rg.vi(4):rg.vf(4))+eFrc(end,:);
            edgeN{3}{feId(edc)}=edgeN{3}{feId(edc)}+eFrc;        
            edgeN{2}{feId(edc)}=edgeN{2}{feId(edc)}+mcp.psi*gmp.dt*eFrc;
%             eFrcAll{feId(edc)}=eFrcAll{feId(edc)}+eFrc;
%             nFrcAll(fac)=nFrcAll(fac)+nFrc;
            faceN{1}(fac,rg.fi(8):rg.ff(8))=...
                faceN{1}(fac,rg.fi(8):rg.ff(8))+nFrc;
            faceN{1}(fac,rg.fi(10):rg.ff(10))=...
                faceN{1}(fac,rg.fi(10):rg.ff(10))+ncMmt;

%             end
        end                 
    end

    % thermal fluctuation on nucleus
    faceN{1}(:,rg.fi(8):rg.ff(8))=faceN{1}(:,rg.fi(8):rg.ff(8))+...
        mcp.mun/mcp.taun*sqrt(gmp.dt).*normrnd(0,1,gmp.nFa,2);
end
    
% nucleus restoring force
cenVec=zeros(gmp.nFa,2);
for fac=1:gmp.nFa
    cenInf=[faceN{1}(fac,rg.fi(6):rg.ff(6));...
        faceN{1}(fac,rg.fi(7):rg.ff(7))];
    cenInf=ns_crdLocal(cenInf,gmp.bs,gmp.sstn);
    cenVec(fac,:)=mcp.nrst*(cenInf(1,:)-cenInf(2,:))/gmp.lsc;
end

faceN{1}(:,rg.fi(8):rg.ff(8))=faceN{1}(:,rg.fi(8):rg.ff(8))+cenVec;

%% Adjust end vertices and update intermediate vertex information
vrtxN(:,rg.vi(3):rg.vf(3))=vrtxN(:,rg.vi(3):rg.vf(3))+...
    gmp.lsc*gmp.dt*vrtxN(:,rg.vi(4):rg.vf(4));
faceN{1}(:,rg.fi(7):rg.ff(7))=faceN{1}(:,rg.fi(7):rg.ff(7))+...
    gmp.lsc*gmp.dt/mcp.taun*faceN{1}(:,rg.fi(8):rg.ff(8));
faceN{1}(:,rg.fi(9):rg.ff(9))=faceN{1}(:,rg.fi(9):rg.ff(9))+...
    gmp.dt/mcp.taua*faceN{1}(:,rg.fi(10):rg.ff(10));

for ii=1:eMx
    eMd=edgeN{2}{ii};
    eMd(1,:)=vrtxN(edgeN{1}(ii,1),rg.vi(3):rg.vf(3));
    eMd(end,:)=vrtxN(edgeN{1}(ii,2),rg.vi(3):rg.vf(3));
    edgeN{2}{ii}=ns_crdLocal(eMd,gmp.bs,gmp.sstn);        
end    

%% Update edge tension.
% [eTn,eTp]=deal(zeros(size(edge{1},1),1));
% for ii=1:eMx
%     [eTn(ii),eTp(ii)]=ns_edgeFixedTension(faceN,edgeN,rg,mcp,ii);
% end
% edgeN{1}(:,rg.ei(3))=edgeN{1}(:,rg.ei(3))...
%     -gmp.dt/mcp.taut*(edgeN{1}(:,rg.ei(3))-eTn)...
%     +mcp.mu/mcp.taut*sqrt(gmp.dt)/2*eTp...
%     .*normrnd(0,1,eMx,1);

edgeN{1}(:,rg.ei(3))=edgeN{1}(:,rg.ei(3))...
    -gmp.dt/mcp.taut*(edgeN{1}(:,rg.ei(3))-1)...
    +mcp.mu*sqrt(2/mcp.taut)*sqrt(gmp.dt)...
    .*normrnd(0,1,eMx,1);
%     
%% Check refinement and adjust number of intermediate vertices
for edc=1:eMx
    emd=edgeN{2}{edc};
    eln=ns_edgeLen(emd);             
    if sum(eln)<gmp.edpc*(edgeN{1}(edc,rg.ei(4))-1/3) || ...
            sum(eln)>gmp.edpc*(edgeN{1}(edc,rg.ei(4))+4/3) || ...
            (max(eln)/min(eln)>2)                
        [edgeN{2}{edc},edgeN{1}(edc,rg.ei(4))]=...
            ns_edgeMidVrtxAverage(emd,gmp.edpc);
        edgeN{4}{edc}=ns_edgeLen(edgeN{2}{edc});
    else
        edgeN{4}{edc}=eln;
    end
    edgeN{1}(edc,rg.ei(5):rg.ef(5))=ns_edgeRange(edgeN,edc);
end

%% calculate area, perimeter, and number of neighbors and center.
for fac=1:gmp.nFa
    fEd=faceN{3}{fac};

    [faceN{1}(fac,rg.fi(2)),faceN{1}(fac,rg.fi(6):rg.ff(6))]...
        =ns_faceArea(fEd,edgeN{2},gmp.bs,gmp.sstn);
    faceN{1}(fac,rg.fi(3))=0;
    for jj=1:size(fEd,2)
        faceN{1}(fac,rg.fi(3))=faceN{1}(fac,rg.fi(3))+...
            sum(edgeN{4}{abs(fEd(jj))});
    end
    faceN{1}(fac,rg.fi(4))=...
        ns_faceNeiCount(fEd,edgeN{1}(:,rg.ei(2):rg.ef(2)),fac);
    faceN{1}(fac,rg.fi(5))=size(fEd,2);
end

vrtxN=ns_crdDomain(vrtxN,rg,gmp);

% toc
end