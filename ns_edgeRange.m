%% Compute edge range.
function erng=ns_edgeRange(edge,eId)

emd=edge{2}{eId};
erng=[min(emd(:,1)),max(emd(:,1)),min(emd(:,2)),max(emd(:,2))];

end