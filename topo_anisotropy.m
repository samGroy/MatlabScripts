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
%Step 5: find c difference beteen source point and all distl points
%Step 6: average difference for each angle, then calculate anisotropy by
%angle
maxradius=1000;
ani=(0:5:175)';
anisotropy=zeros(length(dat(:,1)),2);

for i=1:length(dat(:,3))
    if dat(i,1)<min(dat(:,1))+maxradius || dat(i,1)>max(dat(:,1))-maxradius || dat(i,2)<min(dat(:,2))+maxradius || dat(i,2)>max(dat(:,2))-maxradius
        continue
%     elseif dat(i,3)==0
%         anisotropy(i,1)=-1;
%         anisotropy(i,2)=-1;
%         continue
    else
        dx=dat(i,1)-dat(:,1);
        dy=dat(i,2)-dat(:,2);
        lengthraw=(dx.^2+dy.^2).^.5;
        angleraw=rad2deg(atan((dat(:,2)-dat(i,2))./(dat(:,1)-dat(i,1))));
        angle=angleraw(lengthraw<=maxradius);
        angle(angle<0)=angle(angle<0)+180;
        angle=5*round(angle/5);
        angle(angle==180)=0;
        L=lengthraw(lengthraw<=maxradius);
        c=dat(lengthraw<=maxradius,3);
        ds=abs(c-dat(i,3));
        for j=1:36
            ani(j,2)=mean(ds(angle==ani(j,1)));
            ani((ani(:,2)==0),2)=1e-15;
        end
        anisotropy(i,1)=1-(min(ani(:,2))/max(ani(:,2)));
        if anisotropy(i,1)==0
            anisotropy(i,1:2)=-1;
        else
            anisotropy(i,2)=mean(ani(ani(:,2)==min(ani(:,2)),1));
        end
    end
end