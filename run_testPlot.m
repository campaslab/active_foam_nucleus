ns_plot(edge,face,rg,gmp,mcp);

%% Compute nucleus forces for each edge
eqDis=sqrt(mcp.nvf*mcp.psi/pi);

nFrcAll=zeros(gmp.nFa,2);
eFrcAll=cell(size(edge{1},1),1);
for edc=1:3*gmp.nFa
    eFrcAll{edc}=zeros(size(edge{2}{edc}));
end


for fac=1:gmp.nFa
    feId=abs(face{3}{fac});
    fcen=face{1}(fac,rg.fi(7):rg.ff(7));
    for edc=1:size(feId,2)
        eMd=edge{2}{feId(edc)};
        ecDis=[fcen;eMd];
        ecDis=ns_crdLocal(ecDis,gmp.bs,gmp.sstn);
        fcen=ecDis(1,:);
        eMd=ecDis(2:end,:);
        ecDis=sqrt(sum((ecDis(2:end,:)-ecDis(1,:)).^2,2));
        
        
        if min(ecDis)<eqDis
            [eFrc,nFrc]=ns_nucleusForce(gmp,mcp,fcen,eMd,ecDis);
%             vrtxN(edgeN{1}(feId(edc),1),rg.vi(4):rg.vf(4))=...
%                 vrtxN(edgeN{1}(feId(edc),1),rg.vi(4):rg.vf(4))+eFrc(1,:);
%             vrtxN(edgeN{1}(feId(edc),2),rg.vi(4):rg.vf(4))=...
%                 vrtxN(edgeN{1}(feId(edc),2),rg.vi(4):rg.vf(4))+eFrc(end,:);
%             edgeN{3}{feId(edc)}=edgeN{3}{feId(edc)}+eFrc;        
%             edgeN{2}{feId(edc)}=edgeN{2}{feId(edc)}+mcp.psi*gmp.dt*eFrc;
            eFrcAll{feId(edc)}=eFrcAll{feId(edc)}+eFrc;
            nFrcAll(fac,:)=nFrcAll(fac,:)+nFrc;
%             faceN{1}(fac,rg.fi(8):rg.ff(8))=...
%                 faceN{1}(fac,rg.fi(8):rg.ff(8))+nFrc;
        end
    end                 
end


for edc=1:3*gmp.nFa
    for ptc=1:size(edge{2}{edc},1)
        plot([edge{2}{edc}(ptc,1),edge{2}{edc}(ptc,1)+5*eFrcAll{edc}(ptc,1)],...
            [edge{2}{edc}(ptc,2),edge{2}{edc}(ptc,2)+5*eFrcAll{edc}(ptc,2)],'r',...
            'LineWidth',2);
        hold on;
    end
end

for fac=1:gmp.nFa
    plot([face{1}(fac,rg.fi(7)),face{1}(fac,rg.fi(7))+5*nFrcAll(fac,1)],...
        [face{1}(fac,rg.ff(7)),face{1}(fac,rg.ff(7))+5*nFrcAll(fac,2)],'b',...
        'LineWidth',2);
    hold on;
end