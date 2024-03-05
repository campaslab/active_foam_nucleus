%% Correct vertex coordinates. 
% If vertex crosses the boundary, correct it. 

function vrF=ns_crdLocal(vrL,bs,sstn)

% final vertex coordinats
vrF=vrL;

% correct vertex coordinate if it crosses boundary. 
for j=2:size(vrL,1)
    % adjust y coordinates
    if vrF(j,2)-vrF(j-1,2)>bs/2
        vrF(j,2)=vrF(j,2)-bs;
        vrF(j,1)=vrF(j,1)-bs*sstn;
    end
    if vrF(j,2)-vrF(j-1,2)<-bs/2
        vrF(j,2)=vrF(j,2)+bs;
        vrF(j,1)=vrF(j,1)+bs*sstn;
    end
    
    % adjust x coordinates    
    if vrF(j,1)-vrF(j-1,1)>bs/2
        vrF(j,1)=vrF(j,1)-bs;
    end
    if vrF(j,1)-vrF(j-1,1)<-bs/2
        vrF(j,1)=vrF(j,1)+bs;
    end
end

end