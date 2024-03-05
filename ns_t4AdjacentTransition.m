%% Resolve overlapping vertex
function [vrtxN,edgeN,faceN]=...
    ns_t4AdjacentTransition(vrtx,edge,face,rg,gmp,mcp,vId)
    
%% Define output variables.
[vrtxN,edgeN,faceN]=deal(vrtx,edge,face);

eId=vrtxN(vId,rg.vi(1):rg.vf(1));
eId=eId(eId~=0);

if size(eId,2)>1 
    %% Check whether there is an overlapping
    % list all possible pair of edges that are not adjacent each other
    eId=[eId,eId(1)];
    for ii=1:size(eId,2)-1
        [eCrs,ePt,cr1,cr2]=...
            ns_edgeCrossAdjacent(edgeN,rg,[eId(ii),eId(ii+1)],gmp);
        cr=[cr1,cr2];
        chk=(edgeN{1}(eId(ii),rg.ei(1):rg.ef(1))==...
            edgeN{1}(eId(ii+1),rg.ei(1):rg.ef(1)));
        chk=prod(chk);


        if (eCrs==1 && chk~=1 && ePt(1)~=Inf && ePt(2)~=Inf)
            vrtxN(vId,rg.vi(3):rg.vf(3))=ePt;

            for jj=ii:ii+1
                emd=edgeN{2}{eId(jj)};
                if edgeN{1}(eId(jj),rg.ei(1))==vId
                    emdf=[ePt;emd(cr(jj-ii+1)+1:end,:)];
                else
                    emdf=[emd(1:cr(jj-ii+1),:);ePt];
                end

                emdf=ns_crdLocal(emdf,gmp.bs,gmp.sstn);
                [edgeN{2}{eId(jj)},edgeN{1}(eId(jj),rg.ei(4))]=...
                    ns_edgeMidVrtxAverage(emdf,gmp.edpc);
                edgeN{4}{eId(jj)}=ns_edgeLen(edgeN{2}{eId(jj)});
                edgeN{1}(eId(jj),rg.ei(5):rg.ef(5))=...
                    ns_edgeRange(edgeN,eId(jj));

            end

            fedr=vrtxN(vId,rg.vi(1):rg.vf(1));
            fedr=setdiff(fedr,[eId(ii),eId(ii+1)]);
            emd=edgeN{2}{fedr};
            if edgeN{1}(fedr,rg.ei(1))==vId
                emd(1,:)=vrtxN(vId,rg.vi(3):rg.vf(3));
            else
                emd(end,:)=vrtxN(vId,rg.vi(3):rg.vf(3));
            end
            emd=ns_crdLocal(emd,gmp.bs,gmp.sstn);
            [edgeN{2}{fedr},edgeN{1}(fedr,rg.ei(4))]=...
                    ns_edgeMidVrtxAverage(emd,gmp.edpc);
            edgeN{4}{fedr}=ns_edgeLen(edgeN{2}{fedr});
            edgeN{1}(fedr,rg.ei(5):rg.ef(5))=ns_edgeRange(edgeN,fedr);

            flId=vrtx(vId,rg.vi(2):rg.vf(2));
            flId=flId(flId~=0);
            for kk=1:size(flId,2)
                fed=faceN{3}{flId(kk)};

                [faceN{1}(flId(kk),rg.fi(2)),faceN{1}(flId(kk),rg.fi(6):rg.ff(6))]...
                    =ns_faceArea(fed,edgeN{2},gmp.bs,gmp.sstn);
                faceN{1}(flId(kk),rg.fi(3))=0;
                for jj=1:size(fed,2)
                    faceN{1}(flId(kk),rg.fi(3))=faceN{1}(flId(kk),rg.fi(3))+...
                        sum(edgeN{4}{abs(fed(jj))});
                end
            end

            elnOld=[sum(edge{4}{eId(ii)}),sum(edge{4}{eId(ii+1)})];
            elnNew=[sum(edgeN{4}{eId(ii)}),sum(edgeN{4}{eId(ii+1)})];

            cri=min(elnNew./elnOld);

            far=min(faceN{1}(flId,rg.fi(2)));
            if cri<0.05 || far<10^-2
                [vrtxN,edgeN,faceN]=deal(vrtx,edge,face);
            end                    
        end
    end  
end

end