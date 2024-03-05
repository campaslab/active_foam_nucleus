%% Check whether two edges are crossing or not.
% Check crossing of two adjacent edges.
function [ecr,ept,crps1,crps2]=ns_edgeCrossAdjacent(edge,rg,eId,gm_p)

emd1=edge{2}{eId(1)};
emd2=edge{2}{eId(2)};

emd=[emd1;emd2];
emd=ns_crdLocal(emd,gm_p.bs,gm_p.sstn);

emd1=emd(1:size(emd1,1),:);
emd2=emd(size(emd1,1)+1:end,:);

ecrCk=zeros(size(emd1,1)-1,size(emd2,1)-1);
ept=cell(size(emd1,1)-1,size(emd2,1)-1);

evr1=edge{1}(eId(1),rg.ei(1):rg.ef(1));
evr2=edge{1}(eId(2),rg.ei(1):rg.ef(1));

if evr1(1)==evr2(1)
    eadj=[1,1];
elseif evr1(1)==evr2(2)
    eadj=[1,size(emd2,1)-1];
elseif evr1(2)==evr2(1)
    eadj=[size(emd1,1)-1,1];
else
    eadj=[size(emd1,1)-1,size(emd2,1)-1];
end

[crps1,crps2]=deal(zeros((size(emd1,1)-1)*(size(emd2,1)-1),1));
psc=1;

for ii=1:size(emd1,1)-1
    for jj=1:size(emd2,1)-1
        if ii~=eadj(1) || jj~=eadj(2)
            [ecrCk(ii,jj),pt]=...
                ns_lineCross([emd1(ii:ii+1,:);emd2(jj:jj+1,:)]);
            if ecrCk(ii,jj)==1
                ept{ii,jj}=pt;
                crps1(psc)=ii;
                crps2(psc)=jj;
            end
        end
        psc=psc+1;
    end
end

ecr=sum(ecrCk(:));
ept=cell2mat(ept(:));
crps1=crps1(crps1~=0);
crps2=crps2(crps2~=0);

end