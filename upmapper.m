function upmap=upmapper(basenm,upmin,upmax,outlet_id)
% Simple function that creates an uplift map with gradient increaseing from
% 0,0 to max(x),max(y). Perhaps later I can think up something more clever
% for oblique uplift.

% Read child prelim file, get x,y coords
cname = basenm;
xyzb = creadxyzb(cname,1);
nodes=length(xyzb);
for i=1:nodes
        xc(i) = xyzb(i,1);
        yc(i) = xyzb(i,2);
end

if strcmp(outlet_id, 'corner')
%define 3 points
p01=[0,0,upmin]; %origin point
p11=[max(xc),0,(upmin+upmax)/2]; %point taken down strike
p21=[max(xc),max(yc),upmax]; %point taken up dip
end

if strcmp(outlet_id, 'edge')
%define 3 points
p01=[0,0,upmin]; %origin point
p11=[max(xc),0,upmin]; %point taken down strike
p21=[0,max(yc),upmax]; %point taken up dip
end

%define pole orthogonal to plane from strike and dip inputs 
vector011=p11-p01;
vector121=p21-p01;
ortho1=cross(vector011,vector121);
upmap=-((ortho1(1)*(xc-p01(1))+ortho1(2)*(yc-p01(2)))/ortho1(3));
upmap=upmap';

%viz
tri=delaunay(xc,yc);
trisurf(tri,xc,yc,upmap); 

save upmap001 upmap /ascii;