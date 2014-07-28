function upmapper3(basenm,step)
%This creats a simple antiformal uplift field for CHILD
xyzb = creadxyzb(basenm,step);
nodes=length(xyzb);
for i=1:nodes
    xc(i) = xyzb(i,1);
    yc(i) = xyzb(i,2);
    zc(i) = xyzb(i,3);
    zci(i)=zc(i);
    bc(i) = xyzb(i,4);
end

center=mean(yc);
% up=.5*(1-sin(90*((yc-center)/center)));
up=.001*sin(deg2rad(90*(yc/(max(yc)/2))));
save upmap001 up /ascii;