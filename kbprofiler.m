% Given kb matrix from a qmovie run and xp,yp from a cstrmprofile_vid run,
% creates a stream profile of kb.

for i=1:length(xp)
    dx=dat(:,1)-xp(i); dx=dx.^2;
    dy=dat(:,2)-yp(i); dy=dy.^2;
    dist=sqrt(dx+dy);
    [val,ind]=min(dist);
    kb(i)=dat(ind,3);
end