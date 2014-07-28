% Topo_anisotropy: code to determine directional dependence of topography
% over different scales based on elevation.

%THIS ONE GRABS SLOPE DATA FROM ALL POINT WITHIN A RADIUS AND AVERAGES
%ALL SLOPE DIFFERENCES BASED ON ANGLE
%get topo data

%procedure goes point by point:
%take transects that pass through point i, by 5 degree intervals

%procdure goes transect by transect, using least squares identify a distal point of specified distance from source point

%measure curvature (2nd spatial derivative) between source point and distal point
%Curvature will tell of the transgressions in clope you get between valley
%floors, vally walls, and ridgetops.

%Requires curvature threshold to find a "max" length value, after some runs
%find what a good threshold value of curvature might be. When threshold is
%met on either side of source point, take coordinates and determine the
%length between the two distal points, and take note of the orientation.
%Also measure distance and orientation of the two closest distal points.

%anisotropy in this case will be 1 - shortest axis / longest axis * 100.

% a measure like this will take an immense amount of time for hight density
% model topography, as every point needs measure, and every point needs
% 36 transects done, and every transect could involve finding 100 or so distal
% points.

%Test at a single scale to see if anisotropy maps work

%Step 1: run qmovie and get x, y, slope (option 4) usually called dat
%Step 2: get list of all points at least 500m away from the border (or some
%other distance)

%Step 3: for every point, get distance and angle to all other points within
%500m radius
%Step 4: split distal points by angle
%Step 5: find slope difference beteen source point and all distl points
%Step 6: average difference for each angle, then calculate anisotropy by
%angle
basenm='regimetest45_3';
steps=5;

% collect stream network data
filesys=[''];
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.area' ];
afid=fopen(filenm,'r');
filenm= [filesys basenm '.nodes' ];
nfid=fopen(filenm,'r');


for i=1:steps
    t = fscanf(netfid,'%f',1);
    nnint = fscanf(netfid,'%d',1);
    net = fscanf(netfid,'%f',[1,nnint]);
    
    t = fscanf(afid,'%f',1);
    nnint = fscanf(afid,'%d',1);
    a = fscanf(afid,'%f',[1,nnint]);
    
    t = fscanf(nfid,'%f',1);
    nn = fscanf(nfid,'%d',1);
    coords = fscanf(nfid,'%f',[4,nn]);
end

coords=coords';
net=net';
a=a';
angle=zeros(nnint,1);
% net=net+1;

for i=1:nnint
    dy=coords(net(i),2)-coords(i,2);
    dx=coords(net(i),1)-coords(i,1);
%     angle(i)=rad2deg(atan((coords(net(i),2)-coords(i,2))/(coords(net(i),1)-coords(i,1))));
    angle(i)=atan(dy/dx);
    
%    if coords(i,1)<coords(net(i),1) && coords(i,2)<coords(net(i),2) % quadrant I
%        angle(i)=rad2deg(pi/2-abs(atan((coords(net(i),2)-coords(i,2))/(coords(net(i),1)-coords(i,1)))));
%    end
%    if coords(i,1)>coords(net(i),1) && coords(i,2)<coords(net(i),2) % quadrant i
%        angle(i)=rad2deg(abs(atan((coords(net(i),2)-coords(i,2))/(coords(net(i),1)-coords(i,1)))));
%    end
%    if coords(i,1)>coords(net(i),1) && coords(i,2)>coords(net(i),2) % quadrant i
%        angle(i)=rad2deg(pi/2+abs(atan((coords(net(i),2)-coords(i,2))/(coords(net(i),1)-coords(i,1)))));
%    end
%    if coords(i,1)<coords(net(i),1) && coords(i,2)>coords(net(i),2) % quadrant IV
%        angle(i)=rad2deg(2*pi/2-abs(atan((coords(net(i),2)-coords(i,2))/(coords(net(i),1)-coords(i,1)))));
%    end

%     plot(coords(i,2), coords(i,1),'.r',coords(net(i),2), coords(net(i),1),'.k'); hold on;
%     if angle(i)<0
%         angle(i)=angle(i)+180;
%     end
end

% angle=5*round(angle/5);
% angle(angle==180)=0;angle,


polar(angle,a,'.r');