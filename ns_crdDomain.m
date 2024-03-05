%% Compute new vertex coordinates that is outside of box
% Box shape is now parallelogram so any vertex outside the boundary should
% be positioned inside the box using periodic boundary condition.

function vrtxN=ns_crdDomain(vrtx,rg,gmp)

vrtxN=vrtx;
for ii=1:size(vrtxN,1)
    % Adjust y coordinate
    if vrtxN(ii,rg.vf(3))>gmp.bs
        vrtxN(ii,rg.vf(3))=vrtxN(ii,rg.vf(3))-gmp.bs;
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))-gmp.bs*gmp.sstn;
    elseif vrtxN(ii,rg.vf(3))<0
        vrtxN(ii,rg.vf(3))=vrtxN(ii,rg.vf(3))+gmp.bs;
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))+gmp.bs*gmp.sstn;
    end

    % Adjust x coordinate
    if vrtxN(ii,rg.vi(3))>gmp.bs+gmp.sstn*vrtxN(ii,rg.vf(3))
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))-gmp.bs;
    elseif vrtxN(ii,rg.vi(3))<gmp.sstn*vrtxN(ii,rg.vf(3))
        vrtxN(ii,rg.vi(3))=vrtxN(ii,rg.vi(3))+gmp.bs;
    end
end

end