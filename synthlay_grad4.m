function [yabba,dabba,doo] = synthlay_grad4( basenm, x01,y01,z01,x02,y02,z02,x03,y03,z03,x04,y04,z04,x05,y05,z05,x06,y06,z06,x07,y07,z07,x08,y08,z08,x09,y09,z09,x010,y010,z010,x011,y011,z011,dipangle1,dipangle2,dipangle3,dipangle4,dipangle5,dipangle6,dipangle7,dipangle8,dipangle9,dipangle10,dipangle11,strikeangle1,strikeangle2,strikeangle3,strikeangle4,strikeangle5,strikeangle6,strikeangle7,strikeangle8,strikeangle9,strikeangle10,strikeangle11,fault_thickness1,fault_thickness2,fault_thickness3,fault_thickness4,fault_thickness5,fault_thickness6,fault_thickness7,fault_thickness8,fault_thickness9,fault_thickness10,fault_thickness11,c1,c2,c3,c4,ts,step,numg,gfrac1,gfrac2,gfrac3,gfrac4)
% FOR ANTECEDENCE
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
% Added input functions for grain size info.

% SGR 5/2013 THOUGHT
% why not indicate number of faults, then run a for loop instead of having
% to define 11 faults each time? Contain fault data in a single matrix.
% This would shorten the code length and add flexibility to the number of
% faults that can be generated.

% SGR 5/2013 ANTECEDENCE
% Option for flat lying layer at surface that supercedes
% all faults, i.e. faults do not surface IF their z0 values are below the
% surface. The overlying flat layer is as thick as the difference between
% surface and z0 height. All faults must have the same z0 or things get
% messy...

% Convert degree input to rad value
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

% Convert dip angle to slope value, needed for plane calculation
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

% Calculate vertical thickness of fault plane
w_mult=5; %multiplier determines total length of fault core and damage zone, input 'fault_thickness*' is fault core width
fthick1=abs((fault_thickness1*w_mult)/cos(dipangle1)); %vertical thickness of fault plane
fthick2=abs((fault_thickness2*w_mult)/cos(dipangle2));
fthick3=abs((fault_thickness3*w_mult)/cos(dipangle3));
fthick4=abs((fault_thickness4*w_mult)/cos(dipangle4));
fthick5=abs((fault_thickness5*w_mult)/cos(dipangle5));
fthick6=abs((fault_thickness6*w_mult)/cos(dipangle6));
fthick7=abs((fault_thickness7*w_mult)/cos(dipangle7));
fthick8=abs((fault_thickness8*w_mult)/cos(dipangle8));
fthick9=abs((fault_thickness9*w_mult)/cos(dipangle9));
fthick10=abs((fault_thickness10*w_mult)/cos(dipangle10));
fthick11=abs((fault_thickness11*w_mult)/cos(dipangle11));

% XYZ coordinates for origin point
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
int3=-100000; % the bottom of the domain

% For plane calc, point taken down strike from origin
p11=[x01+cos(strikeangle1),y01+sin(strikeangle1),z01];
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

%For plane calc, point taken down dip from origin
p21=[x01+cos(strikeangle1+(pi/2)),y01+sin(strikeangle1+(pi/2)),slope1+z01];
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

%for plane calc, define pole orthogonal to plane from strike and dip inputs
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

% maximum number of possible layers including 11 faults and the spaces
% between
maxnumlay=67;
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
% planes a-f are the boundaries for each sub layer in the fault. Each fault
% has 5 layers with 3 different cohesion values, designated further below.
% This is done for 11 faults.
intnode=length(xc);
zlt1a=((ortho1(1)*(xc-p01(1))+ortho1(2)*(yc-p01(2)))/ortho1(3))+z01;
zlt1b=zlt1a+fthick1;
zlt1c=zlt1a+fthick1/5;
zlt1d=zlt1a+4*fthick1/5;
zlt1e=zlt1a+2*fthick1/5;
zlt1f=zlt1a+3*fthick1/5;

zlt2a=((ortho2(1)*(xc-p02(1))+ortho2(2)*(yc-p02(2)))/ortho2(3))+z02;
zlt2b=zlt2a+fthick2;
zlt2c=zlt2a+fthick2/5;
zlt2d=zlt2a+4*fthick2/5;
zlt2e=zlt2a+2*fthick2/5;
zlt2f=zlt2a+3*fthick2/5;

zlt3a=((ortho3(1)*(xc-p03(1))+ortho3(2)*(yc-p03(2)))/ortho3(3))+z03;
zlt3b=zlt3a+fthick3;
zlt3c=zlt3a+fthick3/5;
zlt3d=zlt3a+4*fthick3/5;
zlt3e=zlt3a+2*fthick3/5;
zlt3f=zlt3a+3*fthick3/5;

zlt4a=((ortho4(1)*(xc-p04(1))+ortho4(2)*(yc-p04(2)))/ortho4(3))+z04;
zlt4b=zlt4a+fthick4;
zlt4c=zlt4a+fthick4/5;
zlt4d=zlt4a+4*fthick4/5;
zlt4e=zlt4a+2*fthick4/5;
zlt4f=zlt4a+3*fthick4/5;

zlt5a=((ortho5(1)*(xc-p05(1))+ortho5(2)*(yc-p05(2)))/ortho5(3))+z05;
zlt5b=zlt5a+fthick5;
zlt5c=zlt5a+fthick5/5;
zlt5d=zlt5a+4*fthick5/5;
zlt5e=zlt5a+2*fthick5/5;
zlt5f=zlt5a+3*fthick5/5;

zlt6a=((ortho6(1)*(xc-p06(1))+ortho6(2)*(yc-p06(2)))/ortho6(3))+z06;
zlt6b=zlt6a+fthick6;
zlt6c=zlt6a+fthick6/5;
zlt6d=zlt6a+4*fthick6/5;
zlt6e=zlt6a+2*fthick6/5;
zlt6f=zlt6a+3*fthick6/5;

zlt7a=((ortho7(1)*(xc-p07(1))+ortho7(2)*(yc-p07(2)))/ortho7(3))+z07;
zlt7b=zlt7a+fthick7;
zlt7c=zlt7a+fthick7/5;
zlt7d=zlt7a+4*fthick7/5;
zlt7e=zlt7a+2*fthick7/5;
zlt7f=zlt7a+3*fthick7/5;

zlt8a=((ortho8(1)*(xc-p08(1))+ortho8(2)*(yc-p08(2)))/ortho8(3))+z08;
zlt8b=zlt8a+fthick8;
zlt8c=zlt8a+fthick8/5;
zlt8d=zlt8a+4*fthick8/5;
zlt8e=zlt8a+2*fthick8/5;
zlt8f=zlt8a+3*fthick8/5;

zlt9a=((ortho9(1)*(xc-p09(1))+ortho9(2)*(yc-p09(2)))/ortho9(3))+z09;
zlt9b=zlt9a+fthick9;
zlt9c=zlt9a+fthick9/5;
zlt9d=zlt9a+4*fthick9/5;
zlt9e=zlt9a+2*fthick9/5;
zlt9f=zlt9a+3*fthick9/5;

zlt10a=((ortho10(1)*(xc-p010(1))+ortho10(2)*(yc-p010(2)))/ortho10(3))+z010;
zlt10b=zlt10a+fthick10;
zlt10c=zlt10a+fthick10/5;
zlt10d=zlt10a+4*fthick10/5;
zlt10e=zlt10a+2*fthick10/5;
zlt10f=zlt10a+3*fthick10/5;

zlt11a=((ortho11(1)*(xc-p011(1))+ortho11(2)*(yc-p011(2)))/ortho11(3))+z011;
zlt11b=zlt11a+fthick11;
zlt11c=zlt11a+fthick11/5;
zlt11d=zlt11a+4*fthick11/5;
zlt11e=zlt11a+2*fthick11/5;
zlt11f=zlt11a+3*fthick11/5;

% The base of the domain must also be defined with a plane: zltn sits
% horizontally at depth int3.
zltn=repmat(int3,length(xc),1);

% Introduce emty matrices that will hold data for sublayer thickness,
% cohesion, number of layers at each node, depth to the top of each layer,
% and a number indicating each sublayer of the fault, respectively. cdat
% allows for fault cores to supercede other layers with greater cohesion,
% therefore at any given point in 3D, if weak and strong layers overprint,
% the weak layer prevails. cdat also keeps track of each sublayer in each
% fault, so as not to lose track of the sublayer order, especially when
% several faults overlap at a point.
thickness=zeros(maxnumlay, length(xc));
C=nan(maxnumlay,length(xc));
nlayers=zeros(1,length(xc));
depthdat=nan(maxnumlay+1,length(xc));
cdat=depthdat;

%omit any layer sections that would exist above the surface elevation, this
%terminates the plane at the model surface, places limits on the domain,
%otherwise a layer error would occur. Do this for all planes; 11x6=66.
%Change from version 3 to 4: fault sublayers are clipped by the z input for
%the origin point, z must be less than 0 for this version or an error will
%occur.
for i=1:length(xc);
    if zltn(i)>zci(i), zltn(i)=zci(i); end
    if zlt1a(i)>z01, zlt1a(i)=z01; end
    if zlt1b(i)>z01, zlt1b(i)=z01; end
    if zlt2a(i)>z02, zlt2a(i)=z02; end
    if zlt2b(i)>z02, zlt2b(i)=z02; end
    if zlt3a(i)>z03, zlt3a(i)=z03; end
    if zlt3b(i)>z03, zlt3b(i)=z03; end
    if zlt4a(i)>z04, zlt4a(i)=z04; end
    if zlt4b(i)>z04, zlt4b(i)=z04; end
    if zlt5a(i)>z05, zlt5a(i)=z05; end
    if zlt5b(i)>z05, zlt5b(i)=z05; end
    if zlt6a(i)>z06, zlt6a(i)=z06; end
    if zlt6b(i)>z06, zlt6b(i)=z06; end
    if zlt7a(i)>z07, zlt7a(i)=z07; end
    if zlt7b(i)>z07, zlt7b(i)=z07; end
    if zlt8a(i)>z08, zlt8a(i)=z08; end
    if zlt8b(i)>z08, zlt8b(i)=z08; end
    if zlt9a(i)>z09, zlt9a(i)=z09; end
    if zlt9b(i)>z09, zlt9b(i)=z09; end
    if zlt10a(i)>z010, zlt10a(i)=z010; end
    if zlt10b(i)>z010, zlt10b(i)=z010; end
    if zlt11a(i)>z011, zlt11a(i)=z011; end
    if zlt11b(i)>z011, zlt11b(i)=z011; end
    
    if zlt1c(i)>z01, zlt1c(i)=z01; end
    if zlt1d(i)>z01, zlt1d(i)=z01; end
    if zlt2c(i)>z02, zlt2c(i)=z02; end
    if zlt2d(i)>z02, zlt2d(i)=z02; end
    if zlt3c(i)>z03, zlt3c(i)=z03; end
    if zlt3d(i)>z03, zlt3d(i)=z03; end
    if zlt4c(i)>z04, zlt4c(i)=z04; end
    if zlt4d(i)>z04, zlt4d(i)=z04; end
    if zlt5c(i)>z05, zlt5c(i)=z05; end
    if zlt5d(i)>z05, zlt5d(i)=z05; end
    if zlt6c(i)>z06, zlt6c(i)=z06; end
    if zlt6d(i)>z06, zlt6d(i)=z06; end
    if zlt7c(i)>z07, zlt7c(i)=z07; end
    if zlt7d(i)>z07, zlt7d(i)=z07; end
    if zlt8c(i)>z08, zlt8c(i)=z08; end
    if zlt8d(i)>z08, zlt8d(i)=z08; end
    if zlt9c(i)>z09, zlt9c(i)=z09; end
    if zlt9d(i)>z09, zlt9d(i)=z09; end
    if zlt10c(i)>z010, zlt10c(i)=z010; end
    if zlt10d(i)>z010, zlt10d(i)=z010; end
    if zlt11c(i)>z011, zlt11c(i)=z011; end
    if zlt11d(i)>z011, zlt11d(i)=z011; end
    
    if zlt1e(i)>z01, zlt1e(i)=z01; end
    if zlt1f(i)>z01, zlt1f(i)=z01; end
    if zlt2e(i)>z02, zlt2e(i)=z02; end
    if zlt2f(i)>z02, zlt2f(i)=z02; end
    if zlt3e(i)>z03, zlt3e(i)=z03; end
    if zlt3f(i)>z03, zlt3f(i)=z03; end
    if zlt4e(i)>z04, zlt4e(i)=z04; end
    if zlt4f(i)>z04, zlt4f(i)=z04; end
    if zlt5e(i)>z05, zlt5e(i)=z05; end
    if zlt5f(i)>z05, zlt5f(i)=z05; end
    if zlt6e(i)>z06, zlt6e(i)=z06; end
    if zlt6f(i)>z06, zlt6f(i)=z06; end
    if zlt7e(i)>z07, zlt7e(i)=z07; end
    if zlt7f(i)>z07, zlt7f(i)=z07; end
    if zlt8e(i)>z08, zlt8e(i)=z08; end
    if zlt8f(i)>z08, zlt8f(i)=z08; end
    if zlt9e(i)>z09, zlt9e(i)=z09; end
    if zlt9f(i)>z09, zlt9f(i)=z09; end
    if zlt10e(i)>z010, zlt10e(i)=z010; end
    if zlt10f(i)>z010, zlt10f(i)=z010; end
    if zlt11e(i)>z011, zlt11e(i)=z011; end
    if zlt11f(i)>z011, zlt11f(i)=z011; end
    
    %     start sussing out layer data
    %     Assign the cdat code to each plane, which serve as the sublayer
    %     boundaries. zlt#b is the top of the fault, zlt#a is the bottom. c-f
    %     lie in the middle. The value 0 indicates the beginning of a new fault
    %     along a column.
    depthdat(1,i)=zci(i); cdat(1,i)=1;
    depthdat(2,i)=zlt1b(i); cdat(2,i)=0;
    depthdat(3,i)=zlt1d(i); cdat(3,i)=2;
    depthdat(4,i)=zlt1f(i); cdat(4,i)=4;
    depthdat(5,i)=zlt1e(i); cdat(5,i)=5;
    depthdat(6,i)=zlt1c(i); cdat(6,i)=3;
    depthdat(7,i)=zlt1a(i); cdat(7,i)=1;
    depthdat(8,i)=zlt2b(i); cdat(8,i)=0;
    depthdat(9,i)=zlt2d(i); cdat(9,i)=2;
    depthdat(10,i)=zlt2f(i); cdat(10,i)=4;
    depthdat(11,i)=zlt2e(i); cdat(11,i)=5;
    depthdat(12,i)=zlt2c(i); cdat(12,i)=3;
    depthdat(13,i)=zlt2a(i); cdat(13,i)=1;
    depthdat(14,i)=zlt3b(i); cdat(14,i)=0;
    depthdat(15,i)=zlt3d(i); cdat(15,i)=2;
    depthdat(16,i)=zlt3f(i); cdat(16,i)=4;
    depthdat(17,i)=zlt3e(i); cdat(17,i)=5;
    depthdat(18,i)=zlt3c(i); cdat(18,i)=3;
    depthdat(19,i)=zlt3a(i); cdat(19,i)=1;
    depthdat(20,i)=zlt4b(i); cdat(20,i)=0;
    depthdat(21,i)=zlt4d(i); cdat(21,i)=2;
    depthdat(22,i)=zlt4f(i); cdat(22,i)=4;
    depthdat(23,i)=zlt4e(i); cdat(23,i)=5;
    depthdat(24,i)=zlt4c(i); cdat(24,i)=3;
    depthdat(25,i)=zlt4a(i); cdat(25,i)=1;
    depthdat(26,i)=zlt5b(i); cdat(26,i)=0;
    depthdat(27,i)=zlt5d(i); cdat(27,i)=2;
    depthdat(28,i)=zlt5f(i); cdat(28,i)=4;
    depthdat(29,i)=zlt5e(i); cdat(29,i)=5;
    depthdat(30,i)=zlt5c(i); cdat(30,i)=3;
    depthdat(31,i)=zlt5a(i); cdat(31,i)=1;
    depthdat(32,i)=zlt6b(i); cdat(32,i)=0;
    depthdat(33,i)=zlt6d(i); cdat(33,i)=2;
    depthdat(34,i)=zlt6f(i); cdat(34,i)=4;
    depthdat(35,i)=zlt6e(i); cdat(35,i)=5;
    depthdat(36,i)=zlt6c(i); cdat(36,i)=3;
    depthdat(37,i)=zlt6a(i); cdat(37,i)=1;
    depthdat(38,i)=zlt7b(i); cdat(38,i)=0;
    depthdat(39,i)=zlt7d(i); cdat(39,i)=2;
    depthdat(40,i)=zlt7f(i); cdat(40,i)=4;
    depthdat(41,i)=zlt7e(i); cdat(41,i)=5;
    depthdat(42,i)=zlt7c(i); cdat(42,i)=3;
    depthdat(43,i)=zlt7a(i); cdat(43,i)=1;
    depthdat(44,i)=zlt8b(i); cdat(44,i)=0;
    depthdat(45,i)=zlt8d(i); cdat(45,i)=2;
    depthdat(46,i)=zlt8f(i); cdat(46,i)=4;
    depthdat(47,i)=zlt8e(i); cdat(47,i)=5;
    depthdat(48,i)=zlt8c(i); cdat(48,i)=3;
    depthdat(49,i)=zlt8a(i); cdat(49,i)=1;
    depthdat(50,i)=zlt9b(i); cdat(50,i)=0;
    depthdat(51,i)=zlt9d(i); cdat(51,i)=2;
    depthdat(52,i)=zlt9f(i); cdat(52,i)=4;
    depthdat(53,i)=zlt9e(i); cdat(53,i)=5;
    depthdat(54,i)=zlt9c(i); cdat(54,i)=3;
    depthdat(55,i)=zlt9a(i); cdat(55,i)=1;
    depthdat(56,i)=zlt10b(i); cdat(56,i)=0;
    depthdat(57,i)=zlt10d(i); cdat(57,i)=2;
    depthdat(58,i)=zlt10f(i); cdat(58,i)=4;
    depthdat(59,i)=zlt10e(i); cdat(59,i)=5;
    depthdat(60,i)=zlt10c(i); cdat(60,i)=3;
    depthdat(61,i)=zlt10a(i); cdat(61,i)=1;
    depthdat(62,i)=zlt11b(i); cdat(62,i)=0;
    depthdat(63,i)=zlt11d(i); cdat(63,i)=2;
    depthdat(64,i)=zlt11f(i); cdat(64,i)=4;
    depthdat(65,i)=zlt11e(i); cdat(65,i)=5;
    depthdat(66,i)=zlt11c(i); cdat(66,i)=3;
    depthdat(67,i)=zlt11a(i); cdat(67,i)=1;
    
    %     cdat for the base of the domain is set to zero.
    depthdat(68,i)=zltn(i); cdat(68,i)=0;
    
    %     Based on fault orientations, the depth to each fault probably won't
    %     fall in numerical order. Order the depthdat matrix in descending
    %     order to keep track of the depth order for each fault and the cdat
    %     code number for each sublayer. This is key to keep the layer values
    %     in order and find which layers overlap one another.
    [Y,I]=sort(depthdat(:,i),'descend');
    depthdat(:,i)=depthdat(I,i);
    cdat(:,i)=cdat(I,i);
    
    %     This algorithm keeps track of all of the sublayers for each fault and
    %     host rock. Thicknesses of each sublayer are measured and put in the
    %     thickness matrix. In the case of overlapping sublayers, the lowest
    %     cohesion sublayer overprints all others, and those with equal cohesion are
    %     consolidated as a single sublayer. Fault cores will therefore never
    %     be overprinted, and host rock will always be overprinted in the
    %     presence of faults.
    n=2; tic=1; tally=0; sally=tally; bally=sally; % tally keeps track of the outer layer of the fault between planes a and b, sally keeps track of the internal layers between c and d, bally keeps track of the fault core between e and f
    while n <= maxnumlay+1 && sum(thickness(:,i))<-int3-100 % for a depth between the top and bottom of the model domain.
        %         if tic==1
        %             thickness(tic,i)=zci(i)-z01;
        %             C(tic,i)=c2;
        %             tic=tic+1;
        %         end
        if cdat(n,i)==0 %top of outer fault layer, bottom of host layer
            if tally==0 && sally==0 && bally==0 % this signifies the end of a host rock layer.
                thickness(tic,i)=depthdat(n-1,i)-depthdat(n,i); % thickness is depth difference between planes
                if tic ==1 % surface layer for antecedence problem
                    C(tic,i)=1E7;
                else
                    C(tic,i)=c1; % assign host rock cohesion value to sublayer
                end
                tic=tic+1; % move on to next layer, transition to a new layer with a new cohesion value and thickness.
            end
            if tally>0 && sally==0 && bally==0 % Indicates the outermost sublayer of a fault. Accumulate the outer sublayer thicknesses.
                thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
                C(tic,i)=c2; % outermost fault sublayer cohesion value
            end
            tally=tally+1; % when tally > 0, the sublayer is at least the outermost layer of a fault.
            n=n+1; continue % move on to the next index in order.
        end
        
        if cdat(n,i)==2 %top of middle fault layer, bottom of outer fault layer
            if sally==0 %terminate the outermost fault sublayer thickness, the lower cohesion middle fault layer supercedes.
                thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
                C(tic,i)=c2; % assign the outermost fault sublayer cohesion value to the terminated layer
                tic=tic+1; % move on to next sublayer
            end
            if sally>0 && bally==0 % accumulate the middle layer thicknesses
                thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
                C(tic,i)=c3; % assign middle fault sublayer cohesion value
            end
            sally=sally+1; % When sally > 0, the sublayer is at least the middle layer of the fault
            n=n+1; continue
        end
        
        if cdat(n,i)==4 %top of core fault layer, bottom of middle fault layer, indicates that fault core is met at this depth, all other layers are superceeded.
            if bally==0 %terminate the middle sublayer fault thickness
                thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
                C(tic,i)=c3;
                tic=tic+1;
            end
            if bally>0 % accumulate the core layer thicknesses
                thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
                C(tic,i)=c4; % assign fault core cohesion value
            end
            bally=bally+1; % when bally > 0, the sublyer is at least the fault core.
            n=n+1; continue
        end
        %         layer bottoms
        if cdat(n,i)==5 %bottom of core layer, top of middle layer
            bally=bally-1; % the bottom side of the fault core
            thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
            C(tic,i)=c4;
            n=n+1;
            if bally==0 % if all overlapping fault cores are terminated at this depth,
                tic=tic+1; %  move to next sublayer thickness
            end
        end
        if cdat(n,i)==3 %bottom of middle layer, top of outer layer
            sally=sally-1; % the bottom side of the middle fault sublayer
            thickness(tic,i)=thickness(tic,i)+depthdat(n-1,i)-depthdat(n,i);
            C(tic,i)=c3;
            n=n+1;
            if sally==0
                tic=tic+1;
            end
        end
        
        if cdat(n,i)==1 %bottom of outer layer, top of host
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

% clip layers for plot, not necessary for calculation, just visualization,
% can be commented out.
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

% plotting the planes, can be commented out as well.
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
% %convert cohesion to erodibility
% Greg's method
% alpha = 100;
% epsilon = 1;
% kb = alpha ./ ( C + epsilon );
% Hanson and Simon, 2001 conversion for Kb using taucrit
kb=0.2*C.^-0.5;

% if any(thickness==0)
%     error('must have layer thickness greater than zero. May have problems with dip angle.');
% end


% Proceeding section edited from Greg's code for making a layer file.
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
        %         if C(layindex(i),j)==c1
        fprintf(lfid,'%.2f %f %.0f\n',thickness(layindex(i),j),kb(layindex(i),j),0); %end
        %         if C(layindex(i),j)==c2 || C(layindex(i),j)==c3
        %             fprintf(lfid,'%.2f %f %.0f\n',thickness(layindex(i),j),kb(layindex(i),j),1); end
        if numg==1
            fprintf(lfid,'%.2f\n',thickness(layindex(i),j));
        end
        if numg==2
            if C(layindex(i),j)==c1;
                fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac1),(thickness(layindex(i),j)*gfrac2)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
            if C(layindex(i),j)==c2 || C(layindex(i),j)==c3
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


