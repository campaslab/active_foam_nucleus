%% function description.
% for a given seed point. Generate Voronoi diagram with periodic boundary
% condition and determine internal vertices and the relation between vertex
% id and face. Also check the boundary crossing of vertices. 

function [vrtx,edge,face,rg]=ns_randomVoronoi(gmp,sdPt)

%% Define output variables
% specify range of data structure for vrtx, edge, and face
rg=struct;
rg.vi=[1,4,7,9];
rg.vf=[3,6,8,10];
rg.ei=[1,3,5,6,7];
rg.ef=[2,4,5,6,10];
rg.fi=[1,2,3,4,5,6,8,10,12,13];
rg.ff=[1,2,3,4,5,7,9,11,12,13];

vrtx=zeros(2*gmp.nFa,rg.vf(end));
edge=cell(4,1);
edge{1}=zeros(3*gmp.nFa,rg.ef(end));
[edge{2},edge{3},edge{4}]=deal(cell(3*gmp.nFa,1));
face=cell(3,1);
face{1}=zeros(gmp.nFa,rg.ff(end));
[face{2},face{3}]=deal(cell(gmp.nFa,1));

%% Voronoi diagram generation    
% copy seed points for periodic boundary condition
spCp=zeros(9*gmp.nFa,2);
spCp(1:gmp.nFa,:)=sdPt;
pb=[1,1;1,0;1,-1;0,1;0,-1;-1,1;-1,0;-1,-1]*gmp.bs;
for ii=1:size(pb,1)
    spCp(gmp.nFa*ii+1:gmp.nFa*(ii+1),:)=...
        sdPt+repmat(pb(ii,:),gmp.nFa,1);
end

% generate Voronoi diagram. Sort vertices for each face in
% counter-clock wise.
[vr,fVrIni]=voronoin(spCp);
for ii=1:gmp.nFa
    fVrIni{ii}=ns_vrtxIdSort(fVrIni{ii},vr,gmp.bs,0);
end
fVrIni=fVrIni.';

%% Identify vertex-face relation for all vertices. 
vFaMat=zeros(9*gmp.nFa,size(vr,1));
for ii=1:9*gmp.nFa
   vFaMat(ii,fVrIni{ii})=1; 
end
ivId=unique(cell2mat(fVrIni(1:gmp.nFa)));

ivFa=zeros(size(ivId,2),3);
for ii=1:size(ivId,2)
    ivFa(ii,:)=find(vFaMat(:,ivId(ii))==1);
    ivFa(ii,:)=mod(ivFa(ii,:),gmp.nFa)...
        +gmp.nFa*(mod(ivFa(ii,:),gmp.nFa)==0);
    ivFa(ii,:)=sort(ivFa(ii,:));
end
vrtx(:,rg.vi(2):rg.vf(2))=unique(ivFa,'rows');

%% Identify face-vertex relation with new vertex ID. 
vIdPr=zeros(size(ivId,2),2);
vIdPr(:,1)=transpose(ivId);
for ii=1:size(ivId,2)
    vIdPr(ii,2)=find(ismember(vrtx(:,rg.vi(2):rg.vf(2)),...
        ivFa(ii,:),'rows')==1);
end

for ii=1:gmp.nFa
    for jj=1:size(fVrIni{ii},2)
        face{2}{ii}(jj)=vIdPr(vIdPr(:,1)==fVrIni{ii}(jj),2);
    end
end

%% Assign vertex coordinates for new vertex ID.
for ii=1:size(vIdPr,1)
    vrtx(vIdPr(ii,2),rg.vi(3):rg.vf(3))=mod(vr(vIdPr(ii,1),:),gmp.bs);
end

%% Identify edge-vertex relation.
eVr=zeros(gmp.nFa*10,2);
edCnt=1;
for ii=1:gmp.nFa
    fvTmp=face{2}{ii};
    fvTmp=fvTmp(fvTmp~=0);
    vrTmp=[fvTmp,fvTmp(1)];
    for jj=1:size(vrTmp,2)-1
        eVr(edCnt,:)=sort([vrTmp(jj),vrTmp(jj+1)]);
        edCnt=edCnt+1;
    end
end
edge{1}(:,rg.ei(1):rg.ef(1))=unique(eVr(1:edCnt-1,:),'rows');        

%% Identify vertex-edge relation.
veCnt=ones(2*gmp.nFa,1);
for ii=1:3*gmp.nFa
    for jj=1:2
        vrtx(edge{1}(ii,jj),veCnt(edge{1}(ii,jj)))=ii;
        veCnt(edge{1}(ii,jj))=veCnt(edge{1}(ii,jj))+1;
    end
end

%% find edge-face relation.
for ii=1:3*gmp.nFa
    edge{1}(ii,rg.ei(2):rg.ef(2))=sort(intersect(vrtx(edge{1}(ii,1),...
        rg.vi(2):rg.vf(2)),vrtx(edge{1}(ii,2),rg.vi(2):rg.vf(2))));
end  

%% Define intermediate point of initial configuration
for ii=1:3*gmp.nFa
    evr=ns_crdLocal(vrtx(edge{1}(ii,rg.ei(1):rg.ef(1)),...
        rg.vi(3):rg.vf(3)),gmp.bs,gmp.sstn);
    [edge{2}{ii},edge{1}(ii,rg.ei(4))]=...
        ns_edgeMidVrtxAverage(evr,gmp.edpc);               
end

%% Identify face-edge relation.
for ii=1:gmp.nFa
    face{3}{ii}=ns_faceEdgeId(edge,face,gmp,rg,ii);
end       

%% Calculate edge length.
for ii=1:3*gmp.nFa 
    edge{4}{ii}=ns_edgeLen(edge{2}{ii});        
end

%% Assign face type.
% Face type 1 is normal cell and face type 2 is empty space. 
face{1}(:,rg.fi(1))=ones(gmp.nFa,1);

%% calculate area, perimeter, and number of neighbors and center.
for ii=1:gmp.nFa
    fEd=face{3}{ii};

    [face{1}(ii,rg.fi(2)),face{1}(ii,rg.fi(6):rg.ff(6))]...
        =ns_faceArea(fEd,edge{2},gmp.bs,gmp.sstn);
    for jj=1:size(fEd,2)
        face{1}(ii,rg.fi(3))=face{1}(ii,rg.fi(3))+...
            sum(edge{4}{abs(fEd(jj))});
    end
    face{1}(ii,rg.fi(4))=...
        ns_faceNeiCount(fEd,edge{1}(:,rg.ei(2):rg.ef(2)),ii);
    face{1}(ii,rg.fi(5))=size(fEd,2);
end

%% calculate edge tension
edge{1}(:,rg.ei(3))=1;

% edge{1}(:,rg.ei(3))=face{1}(edge{1}(:,3),rg.fi(1))+...
%     face{1}(edge{1}(:,4),rg.fi(1))-...
%     mcp.gam*face{1}(edge{1}(:,3),rg.fi(1))...
%     .*face{1}(edge{1}(:,4),rg.fi(1));

%% Calculate edge range
for ii=1:3*gmp.nFa
    edge{1}(ii,rg.ei(5):rg.ef(5))=ns_edgeRange(edge,ii);
end

%% Assign nucleus position as cell centroid.
face{1}(:,rg.fi(7):rg.ff(7))=face{1}(:,rg.fi(6):rg.ff(6));

%% Assign nucleus orientation randomly
face{1}(:,rg.fi(9):rg.ff(9))=rand(gmp.nFa,1)*2*pi;

end