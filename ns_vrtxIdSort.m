%% Sort vertex ID in counter-clockwise direction. 
% This function is used before periodic boundary condition is applied.

function vSrt=ns_vrtxIdSort(vLst,vCrd,bs,sstn)    

%% Assign face vertex coordinates.
vCrd=ns_crdLocal(vCrd(vLst,:),bs,sstn); 

%% Calculate centroid position and angle about x-axis.
cCrd=[mean(vCrd(:,1)),mean(vCrd(:,2))];
angV=atan2(vCrd(:,2)-cCrd(2),vCrd(:,1)-cCrd(1));

%% Sort vertices in the order of angle value.    
[~,idc]=sort(angV);
vSrt=vLst(idc);    

end