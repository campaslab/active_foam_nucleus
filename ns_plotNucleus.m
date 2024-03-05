%% Plot an individual polygon for a given vertex coordinates. 
function ns_plotNucleus(fCen,orn,gmp)

%% color scheme for nucleus.
% 0:light gray  1:light pink          2:light purple
% 3:red         4:green     5:yellow    6:gray      7:blue
% 8:magenta     9:cyan      10:brown    11:purple
col = [1,0.8,0.8];

%% Determine mirror image of polygon.
pb=[1+gmp.sstn,1;1,0;1-gmp.sstn,-1;...
    gmp.sstn,1;0,0;-gmp.sstn,-1;-1+gmp.sstn,1;-1,0;-1-gmp.sstn,-1]*gmp.bs;

vrTmp=cell(size(pb,1),1);
for ii=1:size(pb,1)
    vrTmp{ii}=repmat(fCen+pb(ii,:),40,1)+...
        gmp.major*[cos(orn),sin(orn)].*(cos(0:pi/20:2*pi-pi/20)).'+...
        gmp.minor*[-sin(orn),cos(orn)].*(sin(0:pi/20:2*pi-pi/20)).';
end

%% Plot polygons with patch function.
for ii=1:size(pb,1)
    patch(vrTmp{ii}(:,1),vrTmp{ii}(:,2),col,'FaceAlpha',0.5);
    hold on;
end
    
end