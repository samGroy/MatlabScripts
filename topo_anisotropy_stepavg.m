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

radius=450; radwindow=50; radstep=10; % measure correlation to radius advancing by radstep
angle=((0:5:355)*pi/180)'; % 5 degree separation of spokes.
spoke=angle.*0; % initialize spoke array
chord=zeros(length(angle)/2,2); % initialize chord array
cor=zeros(length(angle),((radius+radwindow)/radstep)-((radius-radwindow)/radstep)); % initialize raw correlogram matrix
mean_norm_cor=cor; % initialize mean normalized correlogram matrix

x=dat(:,1); %topo X coords
y=dat(:,2); %topo Y coords
c=round(dat(:,3));  %topo Z, rounded to integer (important for fine anisotropy in almost flat regions)

% raw point coords. replace with either clickable map or array of points used to measure.
% xx=1982;
% yy=2095;
% cc=7.048;
% xx=1038;
% yy=3989;
% cc=250.1919;
xx=1951;
yy=1500;
cc=250;

% % measure correlation per radius per angle, generates correlogram matrix
for j=1:((radius+radwindow)/radstep)-((radius-radwindow)/radstep) % Step through radii
    rad=(radius-radwindow)+radstep*j; % Advance radius length per iteration
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
end

% % for j=1:radius/radstep % Step through radii
%     rad=radstep*j; % Advance radius length per iteration
%     for i=1:length(angle) % Step through angles, 5 degree intervals
%         xrad=cos(angle(i))*rad; % Find adjacent length (x component)
%         yrad=sin(angle(i))*rad; % Find opposite length (y component)
%         dx=(xx+xrad)-x(:);
%         dy=(yy+yrad)-y(:);
%         dist=sqrt(dx.^2+dy.^2);
%         [val,ind]=min(dist); % Find point closest to the radius value j for each angle i
%         cor(i,j)=((cc-c(ind)).^2)*.5; % Measure correlation and stick it in an angle x radius matrix
% %         cor(i,j)=cc-c(ind);
%     end
% end
% 

% Average correlogram along spoke
% for i=1:length(cor(1,:))
%     if i==1
%         mean_norm_cor(:,i)=cor(:,i);
%     else
%         mean_norm_cor(:,i)=mean(cor(:,1:i),2);
%     end
% end
   
mean_norm_cor=mean(cor,2);

% Normalize correlogram with respect to minimum and maximum correlation
% values. Within reach of radius.
[high,index1]=max(mean_norm_cor);
[low,index2]=min(mean_norm_cor);
% high=((cc-min(c)).^2)*.5; low=0;
mean_norm_cor=1-(mean_norm_cor-low)./(high-low);
% mean_norm_cor=1./(1+mean_norm_cor);
spoke=mean_norm_cor(:,end)*radius;

% Assume symmetry, Chord sums opposing spokes, later averages to assign semimajor and semiminor axes.
chord(:,1)=0:5:175;
for j=1:length(spoke)/2
    chord(j,2)=spoke(j)+spoke(j+36); 
end

% Calculations to find semimajor and semiminor axes, calculate eccentricity
% or aspect ratio (1-sqrt(max.^2/min.^2)
% figure;
% polar(((0:5:355)*pi/180)',spoke); hold on;
[val1,ind]=max(chord(:,2)); %val2=min(chord(:,2)); % using this min option is kinda cheating, min axis is rarely orthogonal to max axis...
if ind <= 18
    val2=chord(ind+18,2);
else
    val2=chord(ind-18,2);
end
aspect_ratio=axes2ecc(val1,val2);

% Convert vals to degrees for ellipse build (recalc back later):
val1_deg=distdim(val1/2,'m','deg');
xx_deg=distdim(xx,'m','deg');
yy_deg=distdim(yy,'m','deg');

% Create and position the ellipse based on aspect ratio, semimajor axis,
% and source point coordinates
[elat,elon]=ellipse1(yy_deg,xx_deg,[val1_deg aspect_ratio],90-chord(ind,1));
e_y=distdim(elat,'deg','m'); % recalc back to meters
e_x=distdim(elon,'deg','m');

% plot the ellipse, make it so it's transposed on the topo map.
% figure;
hold on;
plot3(e_x,e_y,n,'y'); %need to figure out how to convert degrees to meters for this special case. Ellipse points are in degrees, coordinates are in meters.
% axis([min(x) max(x) min(y) max(y)]); % Here in case you aren't plotting over the topo map
