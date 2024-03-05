%% Compute net normal force for a given edge.
function eNf=ns_edgeNormalForce(edge,face,rg,eId,del,eln,lsc)

eFa=edge{1}(eId,rg.ei(2):rg.ef(2));
eNf=0;
for ii=1:2
    if eFa(ii)~=0
        eNf=eNf+del*(1/(face{1}(eFa(ii),rg.fi(2))/lsc^2)-1)...
            *eln/2/lsc*((-1)^(ii));
    end
end

end