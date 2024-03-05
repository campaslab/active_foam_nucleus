%% Insert vertex id between vr1 and vr2
function fvrsN=ns_faceVrtxIdInsert(vr1,vr2,vrin,fvrs)

fvps=sort([find(fvrs==vr1),find(fvrs==vr2)]);
if fvps(1)==1 && fvps(2)==size(fvrs,2)
    fvrsN=[fvrs,vrin];
else
    fvrsN=[fvrs(1:fvps(1)),vrin,fvrs(fvps(2):end)];
end

end