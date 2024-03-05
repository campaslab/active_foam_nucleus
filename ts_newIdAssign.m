%% Reassign geometic ID
function [vrtxN,edgeN,faceN]=ts_newIdAssign(vrtx,edge,face,rg,vVl,eVl)

[vrtxN,edgeN,faceN]=deal(vrtx,edge,face);
        
% face id
for ii=1:size(face{2},1)
    fvr=face{2}{ii};
    fed=face{3}{ii};
    for jj=1:size(fvr,2)
        fvr(jj)=find(vVl==fvr(jj));
        fed(jj)=sign(fed(jj))*find(eVl==abs(fed(jj)));
    end
    faceN{2}{ii}=fvr;
    faceN{3}{ii}=fed;
end

% vertex id
for ii=1:size(vrtx,1)
    ved=vrtxN(ii,rg.vi(1):rg.vf(1));
    for jj=1:3
        if ved(jj)>0
            vrtxN(ii,rg.vi(1)-1+jj)=find(eVl==ved(jj));
        end
    end
end

% edge id
for ii=1:size(edge{1},1)
    for jj=1:2
        edgeN{1}(ii,rg.ei(1)-1+jj)=find(vVl==edgeN{1}(ii,rg.ei(1)-1+jj));
    end
end 

end