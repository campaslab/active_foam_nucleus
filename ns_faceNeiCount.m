%% Compute number of neighbors.
% Extracellular space cell is not counted as number of neighbors. 
function ne=ns_faceNeiCount(fEd,eFa,fId)

% Find near neighbor faces
neFa=eFa(abs(fEd),:);
neFa=unique(neFa);
neFa=neFa(neFa~=fId);
neFa=neFa(neFa~=0);

% Count number of neighbors if they are not not extracellular space
% cell.    
ne=size(neFa,1);

end