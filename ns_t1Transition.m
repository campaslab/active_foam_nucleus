%% Execute a single T1 trnasition. 

function [vrtxN,edgeN,faceN]=...
    ns_t1Transition(vrtx,edge,face,rg,gmp,mcp,t1Id)

%% Define output variables.
[vrtxN,edgeN,faceN]=deal(vrtx,edge,face);  

%% Find vertex ID, edge ID, and face ID of local T1 configuration.
% local edge ID. 
cvId=edgeN{1}(t1Id,rg.ei(1):rg.ef(1));     
eId=[setdiff(vrtxN(cvId(1),rg.vi(1):rg.vf(1)),t1Id),...
    setdiff(vrtxN(cvId(2),rg.vi(1):rg.vf(1)),t1Id)];
ef1=edgeN{1}(eId(1),rg.ei(2):rg.ef(2));
ef1=ef1(ef1~=0);
if size(ef1,2)==1
    eId(1:2)=[eId(2),eId(1)];
    ef1=edgeN{1}(eId(1),rg.ei(2):rg.ef(2));
    ef1=ef1(ef1~=0);
end

ef2=edgeN{1}(eId(4),rg.ei(2):rg.ef(2));
ef2=ef2(ef2~=0);
if size(ef2,2)==1
    eId(3:4)=[eId(4),eId(3)];
    ef2=edgeN{1}(eId(4),rg.ei(2):rg.ef(2));
    ef2=ef2(ef2~=0);
end

if isempty(intersect(ef1,ef2))==0
    eId(3:4)=[eId(4),eId(3)];
end

% local vertex ID. 
vId=zeros(4,1);    
for ii=1:2
    vId(ii)=setdiff(edgeN{1}(eId(ii),rg.ei(1):rg.ef(1)),cvId(1));
end
for ii=3:4
    vId(ii)=setdiff(edgeN{1}(eId(ii),rg.ei(1):rg.ef(1)),cvId(2));
end

% local face ID.
fId=[intersect(edgeN{1}(t1Id,3:4),edgeN{1}(eId(1),3:4)),...
    intersect(edgeN{1}(t1Id,3:4),edgeN{1}(eId(2),3:4)),...
    setdiff(edgeN{1}(eId(1),3:4),edgeN{1}(t1Id,3:4)),...
    setdiff(edgeN{1}(eId(3),3:4),edgeN{1}(t1Id,3:4))];

%% Reassign vertex-edge relation.
if eId(1)~=eId(3)
    vrtxN(cvId(1),rg.vi(1):rg.vf(1))=sort([t1Id,eId(1),eId(3)]);
else
    vrtxN(cvId(1),rg.vi(1):rg.vf(1))=sort([0,t1Id,eId(1)]);
end

if eId(2)~=eId(4)
    vrtxN(cvId(2),rg.vi(1):rg.vf(1))=sort([t1Id,eId(2),eId(4)]);
else
    vrtxN(cvId(2),rg.vi(1):rg.vf(1))=sort([0,t1Id,eId(2)]);
end


%% Reassign vertex-face relation.
vrtxN(cvId(1),rg.vi(2):rg.vf(2))=sort([fId(1),fId(3),fId(4)]);
vrtxN(cvId(2),rg.vi(2):rg.vf(2))=sort([fId(2),fId(3),fId(4)]);                 

%% Reassign edge-vertex relation.
edgeN{1}(eId(2),rg.ei(1):rg.ef(1))=sort([cvId(2),vId(2)]);
edgeN{1}(eId(3),rg.ei(1):rg.ef(1))=sort([cvId(1),vId(3)]);

%% Reassign edge-face relation.
edgeN{1}(t1Id,rg.ei(2):rg.ef(2))=sort([fId(3),fId(4)]);

%% Adjust vertex coordinates for T1 vertices. 
vCrdTmp=ns_crdLocal(vrtxN(cvId,rg.vi(3):rg.vf(3)),gmp.bs,gmp.sstn);
t1CrdTmp=repmat(1/2*sum(vCrdTmp(:)),2)-[vCrdTmp(2,2),vCrdTmp(1,1);...
    vCrdTmp(1,2),vCrdTmp(2,1)];

emd1=edgeN{2}{eId(1)};
if edge{1}(eId(1),1)==vId(1)
    emd1(end,:)=t1CrdTmp(1,:);
    if vId(1)==cvId(2)
        emd1=emd1(2:end,:);
    end            
else
    emd1(1,:)=t1CrdTmp(1,:);
    if vId(1)==cvId(2)
        emd1=emd1(1:end-1,:);
    end
end

emd2=edgeN{2}{eId(2)};
if edge{1}(eId(2),1)==vId(2)
    emd2(end,:)=t1CrdTmp(2,:);
    if vId(2)==cvId(2)
        emd2=emd2(2:end,:);
    end 
else
    emd2(1,:)=t1CrdTmp(2,:);
    if vId(2)==cvId(2)
        emd1=emd1(1:end-1,:);
    end 
end

emdc=ns_crdLocal([emd1;emd2],gmp.bs,gmp.sstn);
emd1=emdc(1:size(emd1,1),:);
emd2=emdc(size(emd1,1)+1:end,:);

ec=0;
for ii=1:size(emd1,1)-1
    for jj=1:size(emd2,1)-1
        [ecCs,~]=ns_lineCross([emd1(ii:ii+1,:);emd2(jj:jj+1,:)]);
        if ecCs==1
            ec=1;
        end
    end
end

if ec==1
    vrtxN(cvId(1),rg.vi(3):rg.vf(3))=t1CrdTmp(2,:);
    vrtxN(cvId(2),rg.vi(3):rg.vf(3))=t1CrdTmp(1,:);
else
    vrtxN(cvId(1),rg.vi(3):rg.vf(3))=t1CrdTmp(1,:);
    vrtxN(cvId(2),rg.vi(3):rg.vf(3))=t1CrdTmp(2,:);
end

%% Update intermediate edge coordinates.
eVr=ns_crdLocal(vrtxN(edgeN{1}(t1Id,rg.ei(1):rg.ef(1)),...
    rg.vi(3):rg.vf(3)),gmp.bs,gmp.sstn);
[edgeN{2}{t1Id},edgeN{1}(t1Id,rg.ei(4))]=...
    ns_edgeMidVrtxAverage(eVr,gmp.edpc);
if eId(1)~=eId(3)
    for ii=[1,3]
        if edge{1}(eId(ii),1)==vId(ii)
            emd=edge{2}{eId(ii)}(1:end-1,:);
            if edgeN{1}(eId(ii),1)==vId(ii)
                emdc=[emd;vrtxN(cvId(1),rg.vi(3):rg.vf(3))];
            else
                emdc=[vrtxN(cvId(1),rg.vi(3):rg.vf(3));flipud(emd)];
            end
        else
            emd=edge{2}{eId(ii)}(2:end,:);
            if edgeN{1}(eId(ii),2)==vId(ii)
                emdc=[vrtxN(cvId(1),rg.vi(3):rg.vf(3));emd];
            else
                emdc=[flipud(emd);vrtxN(cvId(1),rg.vi(3):rg.vf(3))];
            end
        end            
        eMd=ns_crdLocal(emdc,gmp.bs,gmp.sstn);
        [edgeN{2}{eId(ii)},edgeN{1}(eId(ii),rg.ei(4))]...
            =ns_edgeMidVrtxAverage(eMd,gmp.edpc);
    end
else
    emd=edge{2}{eId(1)}(2:end-1,:);
    emdc=[vrtxN(cvId(1),rg.vi(3):rg.vf(3));emd;...
        vrtxN(cvId(1),rg.vi(3):rg.vf(3))];
    eMd=ns_crdLocal(emdc,gmp.bs,gmp.sstn);
    [edgeN{2}{eId(1)},edgeN{1}(eId(1),rg.ei(4))]...
        =ns_edgeMidVrtxAverage(eMd,gmp.edpc);
end

if eId(2)~=eId(4)
    for ii=[2,4]
        if edge{1}(eId(ii),1)==vId(ii)
            emd=edge{2}{eId(ii)}(1:end-1,:);
            if edgeN{1}(eId(ii),1)==vId(ii)
                emdc=[emd;vrtxN(cvId(2),rg.vi(3):rg.vf(3))];
            else
                emdc=[vrtxN(cvId(2),rg.vi(3):rg.vf(3));flipud(emd)];
            end
        else
            emd=edge{2}{eId(ii)}(2:end,:);
            if edgeN{1}(eId(ii),2)==vId(ii)
                emdc=[vrtxN(cvId(2),rg.vi(3):rg.vf(3));emd];
            else
                emdc=[flipud(emd);vrtxN(cvId(2),rg.vi(3):rg.vf(3))];
            end
        end            
        eMd=ns_crdLocal(emdc,gmp.bs,gmp.sstn);
        [edgeN{2}{eId(ii)},edgeN{1}(eId(ii),rg.ei(4))]...
            =ns_edgeMidVrtxAverage(eMd,gmp.edpc);
    end
else
    emd=edge{2}{eId(2)}(2:end-1,:);
    emdc=[vrtxN(cvId(2),rg.vi(3):rg.vf(3));emd;...
        vrtxN(cvId(2),rg.vi(3):rg.vf(3))];
    eMd=ns_crdLocal(emdc,gmp.bs,gmp.sstn);
    [edgeN{2}{eId(2)},edgeN{1}(eId(2),rg.ei(4))]...
        =ns_edgeMidVrtxAverage(eMd,gmp.edpc);
end

%% Calculate new edge length.
edgeN{4}{t1Id}=ns_edgeLen(edgeN{2}{t1Id});
edgeN{1}(t1Id,rg.ei(5):rg.ef(5))=ns_edgeRange(edgeN,t1Id);
for ii=1:4
    edgeN{4}{eId(ii)}=ns_edgeLen(edgeN{2}{eId(ii)});
    edgeN{1}(eId(ii),rg.ei(5):rg.ef(5))=ns_edgeRange(edgeN,eId(ii));
end

%% Reassign face-vertex relation.
% Insert vertex 2 in face 1 and vertex 1 in face 2. 
for ii=[1,2]
    if fId(ii)>0
        fvrInd=faceN{2}{fId(ii)};
        faceN{2}{fId(ii)}=fvrInd(fvrInd~=cvId(mod(ii,2)+1));
    end
end

for ii=[3,4]
    if fId(ii)>0
        fvrInd=faceN{2}{fId(ii)};
        fvrInd=ns_faceVrtxIdInsert(cvId(mod(ii-1,2)+1),vId(ii-1),...
            cvId(mod(ii,2)+1),fvrInd);
        faceN{2}{fId(ii)}=fvrInd;
    end
end


%% Reassign face-edge relation.
for ii=1:4
    if fId(ii)>0
        faceN{3}{fId(ii)}=ns_faceEdgeId(edgeN,faceN,gmp,rg,fId(ii));
    end
end        

%% Calculate new face area, number of neighbors, perimeter.
for ii=1:4
    if fId(ii)>0
        fEd=faceN{3}{fId(ii)};

        [faceN{1}(fId(ii),rg.fi(2)),faceN{1}(fId(ii),rg.fi(6):rg.ff(6))]...
            =ns_faceArea(fEd,edgeN{2},gmp.bs,gmp.sstn);
        if faceN{1}(fId(ii),rg.fi(2))<0
            faceN{2}{fId(ii)}=fliplr(faceN{2}{fId(ii)});
            faceN{3}{fId(ii)}=ns_faceEdgeId(edgeN,faceN,gmp,rg,fId(ii));
            fEd=faceN{3}{fId(ii)};
            [faceN{1}(fId(ii),rg.fi(2)),faceN{1}(fId(ii),rg.fi(6):rg.ff(6))]...
                =ns_faceArea(fEd,edgeN{2},gmp.bs,gmp.sstn);
        end
        faceN{1}(fId(ii),rg.fi(3))=0;
        for jj=1:size(fEd,2)
            faceN{1}(fId(ii),rg.fi(3))=faceN{1}(fId(ii),rg.fi(3))+...
                sum(edgeN{4}{abs(fEd(jj))});
        end
        faceN{1}(fId(ii),rg.fi(4))=ns_faceNeiCount...
            (fEd,edgeN{1}(:,rg.ei(2):rg.ef(2)),fId(ii));
        faceN{1}(fId(ii),rg.fi(5))=size(fEd,2);
    end
end 

%% Compute edge tension for T1 edge.
% eTn=[ns_edgeFixedTension(face,edge,rg,mcp,t1Id),...
%     ns_edgeFixedTension(faceN,edgeN,rg,mcp,t1Id)];
% 
% edgeN{1}(t1Id,rg.ei(3))=edgeN{1}(t1Id,rg.ei(3))-(eTn(1)-eTn(2));

end