% Topo_anisotropy_correlogram: code to determine directional dependence of topography
% over different scales based on elevation.

% Interest is in continuous correlation, I don't want adjacent valleys to
% influence the signal and skew the semiminor axis. Ways to avoid this:
% Average correlation up to the max radius. Easy to do, problem is that
% solution is not devoted to a single wavelength, it is the mean up to a
% specific wavelength. (topo_anisotropy_avg)
% Terminate search at divides, this will require breaking into the .net and
% searching for points with no contributors. Not worth it...
% ignore it, just use the simple step search and take closest value per
% spoke
% Broad step search over wavelength interval per spoke. Makes for more
% robust solution, may be tricky to pull off.
% Terminate search at minimum correlation. Problem is it yields
% scale-dependent results. (this one)

maxradius=1000;
angle=((0:5:355)*pi/180)';
spoke=angle.*0;
chord=zeros(36,2);

x=dat(:,1);
y=dat(:,2);
c=dat(:,3);
% xx=1982;
% yy=2095;
% cc=7.048;
xx=4044;
yy=3800;
cc=67.99;
radius=10; iniradius=radius;
for j=1:maxradius/iniradius
    radius=iniradius*j;
    for i=1:length(angle)
        xrad=cos(angle(i))*radius;
        yrad=sin(angle(i))*radius;
        dx=(xx+xrad)-x(:);
        dy=(yy+yrad)-y(:);
        dist=sqrt(dx.^2+dy.^2);
        [val,ind]=min(dist);
        cor(i,j)=((cc-c(ind)).^2)*.5;
%         cor(i,j)=cc-c(ind);
    end
end

for j=1:length(spoke)
    for i=1:length(cor(1,:))
        if cor(j,i)>=15 || i==length(cor(1,:))
            spoke(j,1)=iniradius*i;
            break
        end
    end
end

chord(:,1)=0:5:175;
for j=1:length(spoke)/2
    chord(j,2)=spoke(j)+spoke(j+36);
end

% figure;
% polar(((0:5:355)*pi/180)',spoke); hold on;
[val1,ind]=max(chord(:,2)); val2=min(chord(:,2));
aspect_ratio=axes2ecc(val1,val2);

% Convert vals to degrees for ellipse build:
val1_deg=distdim(val1/2,'m','deg');
xx_deg=distdim(xx,'m','deg');
yy_deg=distdim(yy,'m','deg');

% Create and position the ellipse based on aspect ratio, semimajor axis,
% and source point coordinates
[elat,elon]=ellipse1(yy_deg,xx_deg,[val1_deg aspect_ratio],90-chord(ind,1));
e_y=distdim(elat,'deg','m');
e_x=distdim(elon,'deg','m');

% plot the ellipse, make it so it's transposed on the topo map.
% figure;
plot3(e_x,e_y,n,'r'); %need to figure out how to convert degrees to meters for this special case. Ellipse points are in degrees, coordinates are in meters.
axis([min(x) max(x) min(y) max(y)]);

% cor2=(cor-min(cor(:)))./(max(cor(:))-min(cor(:)));
% rho=repmat((iniradius:iniradius:maxradius),length(angle),1);
% theta=repmat((angle),1,length(rho(1,:)));
% [X,Y,C]=pol2cart(theta,rho,cor);
% figure;
% contour(X,Y,C,'Fill','on');
% caxis([0 1]);
