% compute edge length for all edge.

function eLnc=ns_edgeLenAll(edge,gmp)

eMx=size(edge{1},1);
eLnc=zeros(eMx,1);
for ii=1:eMx
    eLnc(ii)=sum(edge{4}{ii});
end

end