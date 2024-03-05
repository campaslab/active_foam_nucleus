%% Compute an individual edge length.
% If an edge is straight, calculate the distance between two end points.

function ln=ns_edgeLen(eCrd)

ln=sqrt(sum((eCrd(2:end,:)-eCrd(1:end-1,:)).^2,2));

end