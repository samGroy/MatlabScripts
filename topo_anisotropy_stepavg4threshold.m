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
window=0;
radius=50; radwindow=4500; radstep=10; % measure correlation to radius advancing by radstep
% radius=500; radwindow=2e4; radstep=2; % measure correlation to radius advancing by radstep
angle=((0:5:355)*pi/180)'; % 5 degree separation of spokes.
spoke=angle.*0; % initialize spoke array
% chord=zeros(length(angle)/2,2); % initialize chord array
if window==1
    cor=zeros(length(angle),radius/radstep-((radius-radwindow)/radstep)); % initialize raw correlogram matrix
else
    cor=zeros(length(angle),radius/radstep); % initialize raw correlogram matrix FOR FULL AVERAGING
end
cor_bi=zeros(length(angle)/2,length(cor(1,:)));
mean_norm_cor=cor; % initialize mean normalized correlogram matrix

x=dat(:,1); %topo X coords
y=dat(:,2); %topo Y coords
c=round(dat(:,3));  %topo Z, rounded to integer (important for fine anisotropy in almost flat regions)


% raw point coords. replace with either clickable map or array of points used to measure.
% Bramaputra
% xx=6.86e5;
% yy=3.024e6;
% cc=104;
% Three rivers
% xx=1.088e6;
% yy=2.957e6;
% cc=1284;
% Hinterland three rivers:river2
% xx=9.694e5;
% yy=3.261e6;
% cc=2621;
% Hinterland three rivers:river3
% xx=1.095e6;
% yy=3.244e6;
% cc=2520;
% Ridge south of Y-S and NBGPM
% xx=6.257e5;
% yy=3.227e6;
% cc=4936;
% Ridge W of three rivers, E of Bramaputra
% xx=9.181e5;
% yy=2.908e6;
% cc=1788;
% Little basin between Bramaputra and 3 rivers
% xx=8.623e5;
% yy=2.939e6;
% cc=225;
% Fault north of tsangpo gorge
% xx=7.346e5;
% yy=3.312e6;
% cc=2638;
% Far southeast
% xx=1.241e6;
% yy=2.821e6;
% cc=2633;
% Ridge between Bramaputra and Little basin
% xx=6.977e5;
% yy=2.914e6;
% cc=1774;
% Ridge far Northeast
% xx=1.261e6;
% yy=3.199e6;
% cc=4573;
% Regimetest45_3, big east tributary branch
xx=cursor_info.Position(1,1);
yy=cursor_info.Position(1,2);
cc=cursor_info.Position(1,3);
% Regimetest45_3, peak east
%     xx=3247;
%     yy=1980;
%     cc=216;
% Regimetest45_3, east trib, further east
%     xx=3018;
%     yy=3008;
%     cc=24;

if radius>xx || radius>yy || radius+xx>max(x) || radius+yy>max(y)
    error('Radius exceeds domain limits at prescribed coordinates. Stop thinking outside the box.');
end

if window==1
    % % measure correlation per radius per angle, generates correlogram matrix
    % WINDOW AVRAGING
    for j=1:radius/radstep-((radius-radwindow)/radstep) % Step through radii
        rad=(radius-radwindow)+radstep*j; % Advance radius length per iteration
        for i=1:length(angle) % Step through angles, 5 degree intervals
            xrad=cos(angle(i))*rad; % Find adjacent length (x component)
            yrad=sin(angle(i))*rad; % Find opposite length (y component)
            dx=(xx+xrad)-x(:);
            dy=(yy+yrad)-y(:);
            dist=sqrt(dx.^2+dy.^2);
            [val,ind]=min(dist); % Find point closest to the radius value j for each angle i
            cor(i,j)=((cc-c(ind)).^2)*.5; % Measure correlation and stick it in an angle x radius matrix
            %         cor(i,j)=((cc-c(ind)))*.5; % Measure correlation and stick it in an angle x radius matrix
            %         cor(i,j)=abs(cc-c(ind));
        end
        if j>1
            cor(:,j)=mean(cor(:,1:j),2);
        end
    end
else
    for j=1:radius/radstep % Step through radii, ALL AVERAGING
        rad=radstep*j; % Advance radius length per iteration
        for i=1:length(angle) % Step through angles, 5 degree intervals
            xrad=cos(angle(i))*rad; % Find adjacent length (x component)
            yrad=sin(angle(i))*rad; % Find opposite length (y component)
            dx=(xx+xrad)-x(:);
            dy=(yy+yrad)-y(:);
            dist=sqrt(dx.^2+dy.^2);
            [val,ind]=min(dist); % Find point closest to the radius value j for each angle i
            cor(i,j)=((cc-c(ind)).^2)*.5; % Measure correlation and stick it in an angle x radius matrix
            %         cor(i,j)=cc-c(ind);
        end
        if j>1
            cor(:,j)=mean(cor(:,1:j),2); %Averages step by step up to max point seeking radius.
        end
    end
end

for j=1:length(cor(:,1))/2
    cor_bi(j,:)=(cor(j,:)+cor(j+36,:))./2;
end

% for k=2:length(cor_bi(1,:))
%     cor_bi(:,k)=mean(cor_bi(:,1:k));
% end

[val1,ind]=min(cor_bi(:,end)); % How about just taking the 1-ratio of strongest correlation/orthogonal correlation value?
if ind <= 18
%     dist2=(val1-cor_bi(ind+18,:)).^2;
%     [val2,ind2]=min(dist2); %Don't do this, instead search from closest to furthest, find first time when value gets close
    val2=cor_bi(ind+18,end);
else
%     dist2=(val1-cor_bi(ind-18,:)).^2;
%     [val2,ind2]=min(dist2);
    val2=cor_bi(ind-18,end); %try this to take 1-ratio of strongest cor and orthogonal cor
end
% for i=1:length(dist2)
% [val2,ind2]=min(dist2);

semimajor=radius;
% semiminor=radstep*ind2;
semiminor=radius*(val1/val2);
aspect_ratio_single=1-sqrt(semiminor^2/semimajor^2);
tilt_single=rad2deg(angle(ind));

% Convert vals to degrees for ellipse build (recalc back later):
rad_deg=distdim(radius,'m','deg');
xx_deg=distdim(xx,'m','deg');
yy_deg=distdim(yy,'m','deg');

% Create and position the ellipse based on aspect ratio, semimajor axis,
% and source point coordinates
% [elat,elon]=ellipse1(yy_deg,xx_deg,[rad_deg aspect_ratio_single],90-chord(ind,1));
[elat,elon]=ellipse1(yy_deg,xx_deg,[rad_deg aspect_ratio_single],90-tilt_single);
e_y=distdim(elat,'deg','m'); % recalc back to meters
e_x=distdim(elon,'deg','m');
n=zeros(length(e_y))+(max(c)+10);

% plot the ellipse, make it so it's transposed on the topo map.
% figure;
hold on;
if window==1
    plot3(e_x,e_y,n,'y','LineWidth',2); %need to figure out how to convert degrees to meters for this special case. Ellipse points are in degrees, coordinates are in meters.
else
    plot3(e_x,e_y,n,'r','LineWidth',2); %need to figure out how to convert degrees to meters for this special case. Ellipse points are in degrees, coordinates are in meters.
end
