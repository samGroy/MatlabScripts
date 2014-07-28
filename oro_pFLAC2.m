% function p=oro_p(basenm,t)
% orographic precipitation script
% Using the equation from Roe, 2005 who adapted it from Sawyer, 1956 and
% Smith, 1979
% SGR 10/2013

x=300000;
dx=xf-x;
for x=min(xf):5000:max(xf)
    dx=xf-x;
    i=1;
for y=min(yf):5000:max(yf)
    dy=yf-y;
    difference=sqrt(dx.^2+dy.^2);
    [dud,index]=min(difference);
    dat(i,1)=index;
    dat(i,2)=yf(index);
    dat(i,3)=zf(index);
    y=y-5000;
end
dat(:,4)=220000-dat(:,2);
slope=diff(dat(:,3))./diff(dat(:,4));
for i=2:length(slope)
    if i==length(slope)
        smoothslope(i)=sum(slope(i-1:i))/2;
    else
        smoothslope(i)=sum(slope(i-1:i+1))/3;
    end
end
%Match up max elevation with max slope positions
[dud,eind]=max(dat(:,3));
[dud,sind]=max(smoothslope);
index=abs(sind-eind);
pprofile(:,2)=dat(index:length(dat(:,3)),3);
pprofile(:,1)=dat(1:length(dat(:,2))+1-index,2);

% triangle=delaunay(xyz(:,1),xyz(:,2));
% nntri=length(triangle);
% % FOR loop slope equation for each triangle
% for i=1:nntri
%     vector1=xyz(triangle(i,2),:)-xyz(triangle(i,1),:); % Point 2 minus point 1
%     vector2=xyz(triangle(i,3),:)-xyz(triangle(i,1),:); % Point 3 minus point 1
%     ortho=cross(vector1,vector2);
%     orthoslope(i)=sqrt(ortho(2).^2+ortho(3).^2);
% end
% for i=1:length(xyz)
%     [r,c,v]=find(triangle==i);
%     slope(i)=mean(orthoslope(r));
% end
trig=zeros(nn,1);
flag=trig;
for ii=1:nnint
    if xyz(ii,2)<xyz(net(ii),2), flag(ii,1)=1; end
    if xyz(ii,2)>xyz(net(ii),2), flag(ii,1)=-1; end
   if xyz(ii,1)<xyz(net(ii),1) && xyz(ii,2)<xyz(net(ii),2) % quadrant I
       trig(ii,1)=rad2deg(pi/2-abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1)))));
   end
   if xyz(ii,1)>xyz(net(ii),1) && xyz(ii,2)<xyz(net(ii),2) % quadrant II
       trig(ii,1)=rad2deg(abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1)))));
   end
   if xyz(ii,1)>xyz(net(ii),1) && xyz(ii,2)>xyz(net(ii),2) % quadrant III
       trig(ii,1)=rad2deg(pi/2+abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1)))));
   end
   if xyz(ii,1)<xyz(net(ii),1) && xyz(ii,2)>xyz(net(ii),2) % quadrant IV
       trig(ii,1)=rad2deg(2*pi/2-abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1)))));
   end
   if net(ii)>nnint
       trig(net(ii),1)=trig(ii,1);
       flag(net(ii),1)=flag(ii,1);
   end
%       trig(ii,2)=((xyz(net(ii),2)-xyz(ii,2)).^2+(xyz(net(ii),1)-xyz(ii,1)).^2).^.5; %measures node to node distance
end
% for ii=1:nnint
%     if xyz(ii,2)<xyz(net(ii),2), flag(ii)=1; end
%     if xyz(ii,2)>xyz(net(ii),2), flag(ii)=-1; end
%    if xyz(ii,1)<xyz(net(ii),1) && xyz(ii,2)<xyz(net(ii),2) % quadrant I
%        trig(ii,1)=rad2deg(abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1)))));
%    end
%    if xyz(ii,1)>xyz(net(ii),1) && xyz(ii,2)<xyz(net(ii),2) % quadrant II
%        trig(ii,1)=rad2deg(abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1))))+pi/2);
%    end
%    if xyz(ii,1)>xyz(net(ii),1) && xyz(ii,2)>xyz(net(ii),2) % quadrant III
%        trig(ii,1)=rad2deg(abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1))))+2*pi/2);
%    end
%    if xyz(ii,1)<xyz(net(ii),1) && xyz(ii,2)>xyz(net(ii),2) % quadrant IV
%        trig(ii,1)=rad2deg(abs(atan((xyz(net(ii),2)-xyz(ii,2))/(xyz(net(ii),1)-xyz(ii,1))))+3*pi/2);
%    end
%       trig(ii,2)=((xyz(net(ii),2)-xyz(ii,2)).^2+(xyz(net(ii),1)-xyz(ii,1)).^2).^.5; %measures node to node distance
% end
% 
% for i=1:nnint
%     dx=(xyz(i,1)-xyz(net(i),1));
%     dy=(xyz(i,2)-xyz(net(i),2));
% %     dx=(xyz((net(:)==i),1)-xyz(i,1));
% %     dy=(xyz((net(:)==i),2)-xyz(i,2));
%     theta(i)=rad2deg(atan(dx/dy));
%     if dx < 0 && theta(i) > 0, theta(i) = -theta(i); end
%     if dy > 0 && theta(i) < 0, theta(i) = -theta(i); end
% end

% theta(nnint+1:nn)=0;
% theta=theta';
gradZs=cos(deg2rad(trig(:,1))).*slope;
% gradZs(flag==0)=gradZs(dist==min(dist(flag~=0)));
% flag(flag==0)=flag(dist==min(dist(flag~=0)));
% gradmax=max(gradZs);
% gradmin=min(gradZs);
% gradZnorm=(gradZs)./(gradmax);
for j=1:length(xyz(:,1)) % Algo to "smooth out" the precip signal isotropically
    dx=xyz(:,1)-xyz(j,1);
    dy=xyz(:,2)-xyz(j,2);
    dist=sqrt(dx.^2+dy.^2);
    if flag(j)==0
        gradZs(j)=mode(gradZs(dist==min(dist(flag~=0))));
        flag(j)=mode(flag(dist==min(dist(flag~=0))));
    end
    gradZcollect=gradZs(dist<=smoothfactor);% & flag==flag(j)); %SGR changed to 50k for coupled model, make function of node spacing.
    flagcollect=flag(dist<=smoothfactor);
    peers=length(flagcollect(flagcollect==flag(j)));
    if peers<=1/4*length(gradZcollect) %If the number of flag peers represent less than 1/5 of the total population in the smoothed area
        gZavg(j)=mean(gradZcollect(flagcollect~=0)); %average out those little guys that disagree with their adjacent neighbors.
    else
        gZavg(j)=mean(gradZcollect(flagcollect==flag(j))); %average only the nodes that agree with node j
    end
end
% gZavgnorm=(gZavg-min(gZavg))/(max(gZavg)-min(gZavg));

rho=1.2; %density of water, kg m^-3
u=10; %incoming wind velocity mag, m s^-1
Hm=2000; %km
% Zs=1000;   %import zc data
% slope=20;
% gradZs=1/100; %import slope data (really figure out wind vector and find downwind point and slope to that point...
% gradZs=tan(deg2rad(slope));
w=sin(deg2rad(gradZs))*u;

% For qsat
a=17.67;
b=243.5; %degC
P=1002; %milbar
T=15; %degC, make dependent on elevation
esat=6.122*exp((a*T)/(b+T));
qsat=0.622*(esat/P);
% qsat=8;
mixing=0.011; %same units as qsat, mass/mass, but isn't necessarily the saturated value

% The more complex options are currently commout'd because they generate
% huge difference between adjacent nodes, it's too sensitive to my gradZs
% calc, so for now I'll just use the flag switch. Runtime drastically
% increases when there's a big data scatter.
% Just ran with simplified form on line 133, big homogeneous blocks speed
% up the solution time drastically.
% S(i)=(rho*qsat*w(i)*gradZs(i)*exp(-xyz(i,3)/Hm))*(31556900/1000)%*(24*365)/1000 % Roe, 2005, units m/y
% smithS=gradZs*w*mixing*rho*(31556900/1000) %Originally units mm/sec, converted to m/y Smith, 1979
% for i=1:nn
%     if flag(i) == 1
%         S(i)=gradZnorm(i)*24;
%     else
%         S(i)=gradZnorm(i)*-6;
%     end
% end
% S(flag==1)=.005;
% S(flag~=1)=.001;
% ridgept=xyz(xyz(:,3)==max(xyz(:,3)),2); % Try defining a ridgepoint for now...think of something later...10/2013
% ridgept=mean(ridgept);
% S(xyz(:,2)>=ridgept)=.005; % Values are pretty high for coupled script.
% S(xyz(:,2)<ridgept)=.001;
% trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3),S);
% Smith, 1979 developed a stable model, dependent mostly on change in slope
% Steps to design windward slope application:
% read in xyz, downstream ID, slope data
% get flow orientation x1y1 always downstream node:
% theta=rad2deg(atan((x2-x1)/(y2-y1)))
% theta(y1(:)>=y2(:))=-theta
% gradZs=cos(deg2rad(theta))*slope
% 
% Other trick: T(z)
% Try this with smoothed gradZ:
gZavg=gZavg';

% S=gZavg.*(.01.*(xyz(:,2)./220000));
ygrad=1.*(xyz(:,2)./220000);
S=.005+(ygrad.*gZavg);
%S=.005+factor;%(.01.*(xyz(:,2)./220000)); %SGR set linear dependency on y coordinate, was .0015 for precip2 run, I think 0.005 for precip1.
% S=abs(gZavg(flag==0)).*.001;%(.01.*(xyz(:,2)./220000)); %SGR set linear dependency on y coordinate, was .0015 for precip2 run, I think 0.005 for precip1.
% condensation_map(:,1:2)=xyz(:,1:2);
% condensation_map(:,3)=S;

% Consider the facet function designed by Davis et al '94 for the PRISM
% model. 5 facets: N-S-E-W-Flat, each has a different efect on precip based
% on the wind trajectory. 
% Algo to find ascending y value nodes in x window.
% Get FLAC initial data, nice and orthogonal, find average x and y spacing.
% get ids for nodes that fit between X0 and X1, start measuring their
% precip starting with Y0 and ending with Yn:
% calc windward slope (Z1-Z0/Y1-Y0)
% figure out precip potential, P=P-Rprev
p=S;
tri=delaunay(xyz(:,1),xyz(:,2));
h=trisurf(tri,xyz(:,1),xyz(:,2),xyz(:,3),S);
fclose (netfid);
fclose(trifid);
fclose(sfid);