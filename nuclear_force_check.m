enFrc=cell(3*gmp.nFa,1);
for edc=1:3*gmp.nFa
    enFrc{edc}=zeros(size(edge{2}{edc},1),2);
end

tic
        
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

        [eFrc,nFrc,ncMmt]=ns_nucleusForce(gmp,mcp,fcen,pcen,forn,eMd,ecDis);
        enFrc{feId(edc)}=enFrc{feId(edc)}+eFrc;        
        enFrc{feId(edc)}=enFrc{feId(edc)}+mcp.psi*gmp.dt*eFrc;
    end                 
end
toc

ns_plot(edge,face,rg,gmp);
for edc=1:3*gmp.nFa
    for ptc=1:size(edge{2}{edc},1)
        plot([edge{2}{edc}(ptc,1),edge{2}{edc}(ptc,1)+enFrc{edc}(ptc,1)],[edge{2}{edc}(ptc,2),edge{2}{edc}(ptc,2)+enFrc{edc}(ptc,2)],'r');
        hold on;
    end
end