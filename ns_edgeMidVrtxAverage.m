%% Averaging intermediate vertices
% If intermediate vertices are very close, redistribute them to make almost
% uniform edge length

function [emdn,rfn]=ns_edgeMidVrtxAverage(emd,edpc)

% compute piecewise length
eln=sqrt(sum((emd(2:end,:)-emd(1:end-1,:)).^2,2));

% compute number of intermediate vertices
rfn=floor(sum(eln)/edpc);
eln=eln/sum(eln);

% specify line parameter for piecewise linear function
erg=zeros(size(eln,1)+1,1);
for ii=2:size(erg)
    erg(ii)=sum(eln(1:ii-1));
end

% redistribute intermediate vertices
cnp=(1/(rfn+1))*(1:rfn);
nmd=zeros(size(cnp,2),2);
for ii=1:size(emd,1)-1
    for jj=1:size(nmd,1)
        nmd(jj,:)=nmd(jj,:)+...
            (emd(ii,:)+(cnp(jj)-erg(ii))/(erg(ii+1)-erg(ii))*...
            (emd(ii+1,:)-emd(ii,:)))*(heaviside(erg(ii+1)-cnp(jj))-...
            heaviside(erg(ii)-cnp(jj)));
    end
end    

emdn=[emd(1,:);nmd;emd(end,:)];

end