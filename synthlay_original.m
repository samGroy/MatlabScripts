function [yabba,dabba,doo] = synthlay_original( basenm, x01,y01,z01,x02,y02,z02,x03,y03,z03,x04,y04,z04,x05,y05,z05,x06,y06,z06,x07,y07,z07,x08,y08,z08,x09,y09,z09,x010,y010,z010,x011,y011,z011,dipangle1,dipangle2,dipangle3,dipangle4,dipangle5,dipangle6,dipangle7,dipangle8,dipangle9,dipangle10,dipangle11,strikeangle1,strikeangle2,strikeangle3,strikeangle4,strikeangle5,strikeangle6,strikeangle7,strikeangle8,strikeangle9,strikeangle10,strikeangle11,fault_thickness1,fault_thickness2,fault_thickness3,fault_thickness4,fault_thickness5,fault_thickness6,fault_thickness7,fault_thickness8,fault_thickness9,fault_thickness10,fault_thickness11,c1,c2,ts,step,numg,gfrac1,gfrac2,gfrac3,gfrac4)
%SGR 7/2012
%used to create MC fracture sets, intersecting faults produced by strain
%that can be eroded preferentially by surface processes.
%faults can cross along strike, need to specify in file, origin
%coordinates, dip, and strike for two planes.

% SGR 9/2012 UPDATE
% Layer thickness algorithm made more efficient and capable of handling
% many more faults, front end allows for 6 faults now, but more could be
% added if desired. All faults are uniform single cohesion.

% SGR 10/2012 UPDATE
% Added input functions for gran size info.
% parameters
% basenm='source';
% c1=4E8; %1E7 produces a pretty weak rock with the hanson & simon equation
% c2=1E4;
% numg=1;
% gfrac1=.1;
% gfrac2=.9;
% gfrac3=.9;
% gfrac4=.1;
% the angle slope
% dipangle1=30; % ENTER DIP ANGLE HERE
% dipangle2=-30;
% dipangle3=-50;
% dipangle4=50;
% dipangle5=-80;
% dipangle6=80;
% strikeangle1=0;
% strikeangle2=0;
% strikeangle3=20;
% strikeangle4=-20;
% strikeangle5=50;
% strikeangle6=-50;
dipangle1=degtorad(dipangle1);
strikeangle1=degtorad(strikeangle1);
dipangle2=degtorad(dipangle2);
strikeangle2=degtorad(strikeangle2);
dipangle3=degtorad(dipangle3);
strikeangle3=degtorad(strikeangle3);
dipangle4=degtorad(dipangle4);
strikeangle4=degtorad(strikeangle4);
dipangle5=degtorad(dipangle5);
strikeangle5=degtorad(strikeangle5);
dipangle6=degtorad(dipangle6);
strikeangle6=degtorad(strikeangle6);
dipangle7=degtorad(dipangle7);
strikeangle7=degtorad(strikeangle7);
dipangle8=degtorad(dipangle8);
strikeangle8=degtorad(strikeangle8);
dipangle9=degtorad(dipangle9);
strikeangle9=degtorad(strikeangle9);
dipangle10=degtorad(dipangle10);
strikeangle10=degtorad(strikeangle10);
dipangle11=degtorad(dipangle11);
strikeangle11=degtorad(strikeangle11);

slope1=tan(dipangle1);
slope2=tan(dipangle2);
slope3=tan(dipangle3);
slope4=tan(dipangle4);
slope5=tan(dipangle5);
slope6=tan(dipangle6);
slope7=tan(dipangle7);
slope8=tan(dipangle8);
slope9=tan(dipangle9);
slope10=tan(dipangle10);
slope11=tan(dipangle11);

% fault_thickness1=400;
% fault_thickness2=400;
% fault_thickness3=400;
% fault_thickness4=400;
% fault_thickness5=400;
% fault_thickness6=400;
fthick1=abs(fault_thickness1/cos(dipangle1)); %vertical thickness of fault plane
fthick2=abs(fault_thickness2/cos(dipangle2)); 
fthick3=abs(fault_thickness3/cos(dipangle3)); 
fthick4=abs(fault_thickness4/cos(dipangle4)); 
fthick5=abs(fault_thickness5/cos(dipangle5)); 
fthick6=abs(fault_thickness6/cos(dipangle6)); 
fthick7=abs(fault_thickness7/cos(dipangle7)); 
fthick8=abs(fault_thickness8/cos(dipangle8)); 
fthick9=abs(fault_thickness9/cos(dipangle9)); 
fthick10=abs(fault_thickness10/cos(dipangle10)); 
fthick11=abs(fault_thickness11/cos(dipangle11)); 

%define 3 points
% x01=5000; y01=3000; 
% x02=5000; y02=7000;
% x03=5000; y03=3000; 
% x04=5000; y04=7000;
% x05=5000; y05=3000; 
% x06=5000; y06=7000;
% z01=-0;
% z02=-0;
int3=-3E12;


p01=[x01,y01,z01]; %origin point
p02=[x02,y02,z02];
p03=[x03,y03,z03];
p04=[x04,y04,z04];
p05=[x05,y05,z05];
p06=[x06,y06,z06];
p07=[x07,y07,z07];
p08=[x08,y08,z08];
p09=[x09,y09,z09];
p010=[x010,y010,z010];
p011=[x011,y011,z011];

p11=[x01+cos(strikeangle1),y01+sin(strikeangle1),z01]; %point taken down strike
p12=[x02+cos(strikeangle2),y02+sin(strikeangle2),z02];
p13=[x03+cos(strikeangle3),y03+sin(strikeangle3),z03];
p14=[x04+cos(strikeangle4),y04+sin(strikeangle4),z04];
p15=[x05+cos(strikeangle5),y05+sin(strikeangle5),z05];
p16=[x06+cos(strikeangle6),y06+sin(strikeangle6),z06];
p17=[x07+cos(strikeangle7),y07+sin(strikeangle7),z07];
p18=[x08+cos(strikeangle8),y08+sin(strikeangle8),z08];
p19=[x09+cos(strikeangle9),y09+sin(strikeangle9),z09];
p110=[x010+cos(strikeangle10),y010+sin(strikeangle10),z010];
p111=[x011+cos(strikeangle11),y011+sin(strikeangle11),z011];

p21=[x01+cos(strikeangle1+(pi/2)),y01+sin(strikeangle1+(pi/2)),slope1+z01]; %point taken down dip
p22=[x02+cos(strikeangle2+(pi/2)),y02+sin(strikeangle2+(pi/2)),slope2+z02];
p23=[x03+cos(strikeangle3+(pi/2)),y03+sin(strikeangle3+(pi/2)),slope3+z03]; 
p24=[x04+cos(strikeangle4+(pi/2)),y04+sin(strikeangle4+(pi/2)),slope4+z04];
p25=[x05+cos(strikeangle5+(pi/2)),y05+sin(strikeangle5+(pi/2)),slope5+z05];
p26=[x06+cos(strikeangle6+(pi/2)),y06+sin(strikeangle6+(pi/2)),slope6+z06];
p27=[x07+cos(strikeangle7+(pi/2)),y07+sin(strikeangle7+(pi/2)),slope7+z07];
p28=[x08+cos(strikeangle8+(pi/2)),y08+sin(strikeangle8+(pi/2)),slope8+z08]; 
p29=[x09+cos(strikeangle9+(pi/2)),y09+sin(strikeangle9+(pi/2)),slope9+z09];
p210=[x010+cos(strikeangle10+(pi/2)),y010+sin(strikeangle10+(pi/2)),slope10+z010];
p211=[x011+cos(strikeangle11+(pi/2)),y011+sin(strikeangle11+(pi/2)),slope11+z011];

%define pole orthogonal to plane from strike and dip inputs 
vector011=p11-p01;
vector121=p21-p01;
ortho1=cross(vector011,vector121);
vector012=p12-p02;
vector122=p22-p02;
ortho2=cross(vector012,vector122);
vector013=p13-p03;
vector123=p23-p03;
ortho3=cross(vector013,vector123);
vector014=p14-p04;
vector124=p24-p04;
ortho4=cross(vector014,vector124);
vector015=p15-p05;
vector125=p25-p05;
ortho5=cross(vector015,vector125);
vector016=p16-p06;
vector126=p26-p06;
ortho6=cross(vector016,vector126);
vector017=p17-p07;
vector127=p27-p07;
ortho7=cross(vector017,vector127);
vector018=p18-p08;
vector128=p28-p08;
ortho8=cross(vector018,vector128);
vector019=p19-p09;
vector129=p29-p09;
ortho9=cross(vector019,vector129);
vector0110=p110-p010;
vector1210=p210-p010;
ortho10=cross(vector0110,vector1210);
vector0111=p111-p011;
vector1211=p211-p011;
ortho11=cross(vector0111,vector1211);

% ts=0;

% maximum number of possible layers including 6 faults and the spaces
% between
maxnumlay=23;
% Read child prelim file, omit boundary nodes
cname = basenm;
xyzb = creadxyzb(cname,step);
nodes=length(xyzb);
for i=1:nodes
    if xyzb(i,4)==0
        xc(i) = xyzb(i,1);
        yc(i) = xyzb(i,2);
        zc(i) = xyzb(i,3);
        zci(i)=zc(i);
        bc(i) = xyzb(i,4);
    end
end

% Define the fault planes using the CHILD xy coordinates
intnode=length(xc);
zlt1a=((ortho1(1)*(xc-p01(1))+ortho1(2)*(yc-p01(2)))/ortho1(3))+z01;
zlt1b=zlt1a+fthick1;
zlt2a=((ortho2(1)*(xc-p02(1))+ortho2(2)*(yc-p02(2)))/ortho2(3))+z02;
zlt2b=zlt2a+fthick2;
zlt3a=((ortho3(1)*(xc-p03(1))+ortho3(2)*(yc-p03(2)))/ortho3(3))+z03;
zlt3b=zlt3a+fthick3;
zlt4a=((ortho4(1)*(xc-p04(1))+ortho4(2)*(yc-p04(2)))/ortho4(3))+z04;
zlt4b=zlt4a+fthick4;
zlt5a=((ortho5(1)*(xc-p05(1))+ortho5(2)*(yc-p05(2)))/ortho5(3))+z05;
zlt5b=zlt5a+fthick5;
zlt6a=((ortho6(1)*(xc-p06(1))+ortho6(2)*(yc-p06(2)))/ortho6(3))+z06;
zlt6b=zlt6a+fthick6;
zlt7a=((ortho7(1)*(xc-p07(1))+ortho7(2)*(yc-p07(2)))/ortho7(3))+z07;
zlt7b=zlt7a+fthick7;
zlt8a=((ortho8(1)*(xc-p08(1))+ortho8(2)*(yc-p08(2)))/ortho8(3))+z08;
zlt8b=zlt8a+fthick8;
zlt9a=((ortho9(1)*(xc-p09(1))+ortho9(2)*(yc-p09(2)))/ortho9(3))+z09;
zlt9b=zlt9a+fthick9;
zlt10a=((ortho10(1)*(xc-p010(1))+ortho10(2)*(yc-p010(2)))/ortho10(3))+z010;
zlt10b=zlt10a+fthick10;
zlt11a=((ortho11(1)*(xc-p011(1))+ortho11(2)*(yc-p011(2)))/ortho11(3))+z011;
zlt11b=zlt11a+fthick11;
zltn=repmat(int3,length(xc),1);


thickness=zeros(maxnumlay, length(xc));
C=nan(maxnumlay,length(xc));
nlayers=zeros(1,length(xc));
depthdat=nan(maxnumlay+1,length(xc));
cdat=depthdat;

%omit any layer sections that would exist above the surface elevation
for i=1:length(xc);
    if zlt1a(i)>zci(i)
        zlt1a(i)=zci(i);
    end
    if zlt1b(i)>zci(i)
        zlt1b(i)=zci(i);
    end
    if zlt2a(i)>zci(i)
        zlt2a(i)=zci(i);
    end
    if zlt2b(i)>zci(i)
        zlt2b(i)=zci(i);
    end
    if zltn(i)>zci(i)
        zltn(i)=zci(i);
    end
    if zlt3a(i)>zci(i)
        zlt3a(i)=zci(i);
    end
    if zlt3b(i)>zci(i)
        zlt3b(i)=zci(i);
    end
    if zlt4a(i)>zci(i)
        zlt4a(i)=zci(i);
    end
    if zlt4b(i)>zci(i)
        zlt4b(i)=zci(i);
    end
    if zlt5a(i)>zci(i)
        zlt5a(i)=zci(i);
    end
    if zlt5b(i)>zci(i)
        zlt5b(i)=zci(i);
    end
    if zlt6a(i)>zci(i)
        zlt6a(i)=zci(i);
    end
    if zlt6b(i)>zci(i)
        zlt6b(i)=zci(i);
    end
    if zlt7a(i)>zci(i)
        zlt7a(i)=zci(i);
    end
    if zlt7b(i)>zci(i)
        zlt7b(i)=zci(i);
    end
    if zlt8a(i)>zci(i)
        zlt8a(i)=zci(i);
    end
    if zlt8b(i)>zci(i)
        zlt8b(i)=zci(i);
    end
    if zlt9a(i)>zci(i)
        zlt9a(i)=zci(i);
    end
    if zlt9b(i)>zci(i)
        zlt9b(i)=zci(i);
    end
    if zlt10a(i)>zci(i)
        zlt10a(i)=zci(i);
    end
    if zlt10b(i)>zci(i)
        zlt10b(i)=zci(i);
    end
    if zlt11a(i)>zci(i)
        zlt11a(i)=zci(i);
    end
    if zlt11b(i)>zci(i)
        zlt11b(i)=zci(i);
    end

    %     start sussing out layer data
    depthdat(1,i)=zci(i);
    cdat(1,i)=1;
    depthdat(2,i)=zlt1b(i);
    cdat(2,i)=0;
    depthdat(3,i)=zlt1a(i);
    cdat(3,i)=1;
    depthdat(4,i)=zlt2b(i);
    cdat(4,i)=0;
    depthdat(5,i)=zlt2a(i);
    cdat(5,i)=1;
    depthdat(6,i)=zlt3b(i);
    cdat(6,i)=0;
    depthdat(7,i)=zlt3a(i);
    cdat(7,i)=1;
    depthdat(8,i)=zlt4b(i);
    cdat(8,i)=0;
    depthdat(9,i)=zlt4a(i);
    cdat(9,i)=1;
    depthdat(10,i)=zlt5b(i);
    cdat(10,i)=0;
    depthdat(11,i)=zlt5a(i);
    cdat(11,i)=1;
    depthdat(12,i)=zlt6b(i);
    cdat(12,i)=0;
    depthdat(13,i)=zlt6a(i);
    cdat(13,i)=1;
    depthdat(14,i)=zlt7b(i);
    cdat(14,i)=0;
    depthdat(15,i)=zlt7a(i);
    cdat(15,i)=1;
    depthdat(16,i)=zlt8b(i);
    cdat(16,i)=0;
    depthdat(17,i)=zlt8a(i);
    cdat(17,i)=1;
    depthdat(18,i)=zlt9b(i);
    cdat(18,i)=0;
    depthdat(19,i)=zlt9a(i);
    cdat(19,i)=1;
    depthdat(20,i)=zlt10b(i);
    cdat(20,i)=0;
    depthdat(21,i)=zlt10a(i);
    cdat(21,i)=1;
    depthdat(22,i)=zlt11b(i);
    cdat(22,i)=0;
    depthdat(23,i)=zlt11a(i);
    cdat(23,i)=1;
    depthdat(24,i)=zltn(i);
    cdat(24,i)=0;
    [Y,I]=sort(depthdat(:,i),'descend');
    depthdat(:,i)=depthdat(I,i);
    cdat(:,i)=cdat(I,i);
    
    n=2; tic=1; tally=0;
    while n <= maxnumlay+1
        if cdat(n,i)==0
            if tally==0
                thickness(tic,i)=depthdat(n-1,i)-depthdat(n,i);
                C(tic,i)=c1;
                tic=tic+1;
            end
            if tally>0
                thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
                C(tic,i)=c2;
            end
            tally=tally+1;
            n=n+1; continue
        end
        if cdat(n,i)==1
            tally=tally-1;
            thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
            C(tic,i)=c2;
            n=n+1;
            if tally==0
                tic=tic+1;
            end
        end
    end

end
% round thicknesses to prevent zero layer thickness error
thickness=round(10*thickness)/10;
% plot
for i=1:length(zci)
    if zlt1a(i)==zci(i)
        zlt1a(i)=NaN;
    end
    if zlt1b(i)==zci(i)
        zlt1b(i)=NaN;
    end
    if zlt2a(i)==zci(i)
        zlt2a(i)=NaN;
    end
    if zlt2b(i)==zci(i)
        zlt2b(i)=NaN;
    end
    if zlt3a(i)==zci(i)
        zlt3a(i)=NaN;
    end
    if zlt3b(i)==zci(i)
        zlt3b(i)=NaN;
    end
    if zlt4a(i)==zci(i)
        zlt4a(i)=NaN;
    end
    if zlt4b(i)==zci(i)
        zlt4b(i)=NaN;
    end
    if zlt5a(i)==zci(i)
        zlt5a(i)=NaN;
    end
    if zlt5b(i)==zci(i)
        zlt5b(i)=NaN;
    end
    if zlt6a(i)==zci(i)
        zlt6a(i)=NaN;
    end
    if zlt6b(i)==zci(i)
        zlt6b(i)=NaN;
    end
    if zlt7a(i)==zci(i)
        zlt7a(i)=NaN;
    end
    if zlt7b(i)==zci(i)
        zlt7b(i)=NaN;
    end
    if zlt8a(i)==zci(i)
        zlt8a(i)=NaN;
    end
    if zlt8b(i)==zci(i)
        zlt8b(i)=NaN;
    end
    if zlt9a(i)==zci(i)
        zlt9a(i)=NaN;
    end
    if zlt9b(i)==zci(i)
        zlt9b(i)=NaN;
    end
    if zlt10a(i)==zci(i)
        zlt10a(i)=NaN;
    end
    if zlt10b(i)==zci(i)
        zlt10b(i)=NaN;
    end
    if zlt11a(i)==zci(i)
        zlt11a(i)=NaN;
    end
    if zlt11b(i)==zci(i)
        zlt11b(i)=NaN;
    end

end
tri=delaunay(xc,yc);
trisurf(tri,xc,yc,zlt1a); hold on;
trisurf(tri,xc,yc,zlt1b); hold on;
trisurf(tri,xc,yc,zlt2a); hold on;
trisurf(tri,xc,yc,zlt2b); hold on;
trisurf(tri,xc,yc,zlt3a); hold on;
trisurf(tri,xc,yc,zlt3b); hold on;
trisurf(tri,xc,yc,zlt4a); hold on;
trisurf(tri,xc,yc,zlt4b); hold on;
trisurf(tri,xc,yc,zlt5a); hold on;
trisurf(tri,xc,yc,zlt5b); hold on;
trisurf(tri,xc,yc,zlt6a); hold on;
trisurf(tri,xc,yc,zlt6b); hold on;
trisurf(tri,xc,yc,zlt7a); hold on;
trisurf(tri,xc,yc,zlt7b); hold on;
trisurf(tri,xc,yc,zlt8a); hold on;
trisurf(tri,xc,yc,zlt8b); hold on;
trisurf(tri,xc,yc,zlt9a); hold on;
trisurf(tri,xc,yc,zlt9b); hold on;
trisurf(tri,xc,yc,zlt10a); hold on;
trisurf(tri,xc,yc,zlt10b); hold on;
trisurf(tri,xc,yc,zlt11a); hold on;
trisurf(tri,xc,yc,zlt11b); hold on;

shading flat
% plot3(p01(1),p01(2),p01(3),p11(1),p11(2),p11(3),p21(1),p21(2),p21(3),'or'); hold on; 
% plot3(p02(1),p02(2),p02(3),p12(1),p12(2),p12(3),p22(1),p22(2),p22(3),'ob'); hold off;
axis([0 max(xc) 0 max(yc) -1000 0])
% caxis([-10000,'auto']);
view([-45,45])
%convert cohesion to erodibility
% Greg's method
% alpha = 100;
% epsilon = 1;
% kb = alpha ./ ( C + epsilon );
% Hanson and Simon, 2001 conversion for Kb using taucrit
kb=0.2*C.^-0.5;

% if any(thickness==0)
%     error('must have layer thickness greater than zero. May have problems with dip angle.');
% end

%--now, make the lay file--%
% Create layer file and open it for writing (first make sure it doesn't
% already exist)
if ts==0
    filenm= [ basenm '.lay' ];
else
    %     filenm= [ basenm sprintf('.lay%.0f',step) ];
    filenm= [ basenm '.lay' ];
end
lfid=fopen(filenm,'r');
% if lfid>0
%     fprintf('Layer file "%s" already exists.\n',filenm);
%     error('Try a different name.\n');
% end
lfid = fopen(filenm,'w');
if lfid<1
    error('Unable to create layer file.\n');
end

fprintf('now writing lay file.\n');


% Write header information
fprintf(lfid,' %.2f\n',ts);
fprintf(lfid,'%d\n',length(intnode));

% Loop over nodes to write information
for j=1:length(xc)
    nlayers(j)=sum(thickness(:,j)>0);
    layindex=find(thickness(:,j));
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',nlayers(j));
    
    % For each layer at this node, write information
    for i=1:nlayers(j)
        % Basic properties: creation time, recent activity time, exposure
        % time, thickness, erodibility, and regolith/bedrock flag
        fprintf(lfid,'%.2f %.2f %.2f\n',0,0,0);
        if C(layindex(i),j)==c1, fprintf(lfid,'%.2f %f %.0f\n',thickness(layindex(i),j),kb(layindex(i),j),0); end
        if C(layindex(i),j)==c2, fprintf(lfid,'%.2f %f %.0f\n',thickness(layindex(i),j),kb(layindex(i),j),1); end
        if numg==1
        fprintf(lfid,'%.2f\n',thickness(layindex(i),j));
        end
        if numg==2
            if C(layindex(i),j)==c1;
            fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac1),(thickness(layindex(i),j)*gfrac2)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
        if C(layindex(i),j)==c2
            fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac3),(thickness(layindex(i),j)*gfrac4)); end
        end
        %         % Grain size information
        %         if numg==1
        %             fprintf(lfid,'%.2f\n',laydat(j,i,1));   % this line gives 'old' format: if one size, just write thickness
        %         elseif numg>1
        %             fprintf(lfid,'%.2f ',laydat(j,i,1)-sum(laydat(j,i,5:(3+numg))));  % this line gives 'old' format: if more than one, write thickness of size 1 as total thick minus sum of sizes 2-N
        %             for k=2:numg-1
        %                 fprintf(lfid,'%f ',laydat(j,i,4+(k-1)) );
        %             end
        %             fprintf(lfid,'%f\n',laydat(j,i,4+(numg-1)) );
        %         end
    end
    
end
fprintf('done.\n');
yabba=nlayers;
dabba=thickness;
doo=kb;


