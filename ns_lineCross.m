%% Check whether two line segments intersect or not
function [chk,ptc]=ns_lineCross(vr)

[xx,yy]=deal(vr(:,1),vr(:,2));
det=(xx(2)-xx(1))*(yy(4)-yy(3))-(yy(2)-yy(1))*(xx(4)-xx(3));
if det==0
    chk=0;
    ptc=[inf,inf];
else
    mm=((xx(3)-xx(1))*(yy(2)-yy(1))-(yy(3)-yy(1))*(xx(2)-xx(1)))/det;
    kk=((xx(3)-xx(1))*(yy(4)-yy(3))-(yy(3)-yy(1))*(xx(4)-xx(3)))/det;
    if mm>0 && mm<1 && kk>0 && kk<1
        chk=1;
        ptc=vr(1,:)+kk*(vr(2,:)-vr(1,:));
    else 
        chk=0;
        ptc=[inf,inf];
    end
end
            
end