%% Find vertex coordinates of a given face
% include all intermediate vertices

function [fCrd,feId]=ns_faceVrtxCrd(fEd,eMd,bs,sstn)

fCrd=cell(size(fEd,2),1);
feId=cell(size(fEd,2),1);

for ii=1:size(fEd,2)
    if fEd(ii)>0
        fCrd{ii}=eMd{fEd(ii)}(1:end-1,:);
    else
        fCrd{ii}=flipud(eMd{-fEd(ii)}(2:end,:));
    end
    feId{ii}=repmat(abs(fEd(ii)),size(eMd{abs(fEd(ii))},1)-1,1);
end    
fCrd=cell2mat(fCrd);
fCrd=ns_crdLocal(fCrd,bs,sstn);

feId=cell2mat(feId);
end