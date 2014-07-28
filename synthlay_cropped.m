function [yabba,dabba,doo] = synthlay_grad2( basenm, x01,y01,z01,x02,y02,z02,x03,y03,z03,x04,y04,z04,x05,y05,z05,x06,y06,z06,x07,y07,z07,x08,y08,z08,x09,y09,z09,x010,y010,z010,x011,y011,z011,dip_strike(1,1),dip_strike(2,1),dip_strike(3,1),dip_strike(4,1),dip_strike(5,1),dip_strike(6,1),dip_strike(7,1),dip_strike(8,1),dip_strike(9,1),dip_strike(10,1),dip_strike(11,1),dip_strike(1,2),dip_strike(2,2),dip_strike(3,2),dip_strike(4,2),dip_strike(5,2),dip_strike(6,2),dip_strike(7,2),dip_strike(8,2),dip_strike(9,2),dip_strike(10,2),dip_strike(11,2),fault_thickness1,fault_thickness2,fault_thickness3,fault_thickness4,fault_thickness5,fault_thickness6,fault_thickness7,fault_thickness8,fault_thickness9,fault_thickness10,fault_thickness11,c1,c2,c3,c4,ts,step,numg,gfrac1,gfrac2,gfrac3,gfrac4)
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

dip_strike(1,1)=degtorad(dip_strike(1,1));
dip_strike(1,2)=degtorad(dip_strike(1,2));
dip_strike(2,1)=degtorad(dip_strike(2,1));
dip_strike(2,2)=degtorad(dip_strike(2,2));
dip_strike(3,1)=degtorad(dip_strike(3,1));
dip_strike(3,2)=degtorad(dip_strike(3,2));
dip_strike(4,1)=degtorad(dip_strike(4,1));
dip_strike(4,2)=degtorad(dip_strike(4,2));
dip_strike(5,1)=degtorad(dip_strike(5,1));
dip_strike(5,2)=degtorad(dip_strike(5,2));
dip_strike(6,1)=degtorad(dip_strike(6,1));
dip_strike(6,2)=degtorad(dip_strike(6,2));
dip_strike(7,1)=degtorad(dip_strike(7,1));
dip_strike(7,2)=degtorad(dip_strike(7,2));
dip_strike(8,1)=degtorad(dip_strike(8,1));
dip_strike(8,2)=degtorad(dip_strike(8,2));
dip_strike(9,1)=degtorad(dip_strike(9,1));
dip_strike(9,2)=degtorad(dip_strike(9,2));
dip_strike(10,1)=degtorad(dip_strike(10,1));
dip_strike(10,2)=degtorad(dip_strike(10,2));
dip_strike(11,1)=degtorad(dip_strike(11,1));
dip_strike(11,2)=degtorad(dip_strike(11,2));

slope1=tan(dip_strike(1,1));
slope2=tan(dip_strike(2,1));
slope3=tan(dip_strike(3,1));
slope4=tan(dip_strike(4,1));
slope5=tan(dip_strike(5,1));
slope6=tan(dip_strike(6,1));
slope7=tan(dip_strike(7,1));
slope8=tan(dip_strike(8,1));
slope9=tan(dip_strike(9,1));
slope10=tan(dip_strike(10,1));
slope11=tan(dip_strike(11,1));

w_mult=5; %multiplier determines total length of fault core and damage zone, input 'fault_thickness*' is fault core width
fthick1=abs((fault_thickness1*w_mult)/cos(dip_strike(1,1))); %vertical thickness of fault plane
fthick2=abs((fault_thickness2*w_mult)/cos(dip_strike(2,1))); 
fthick3=abs((fault_thickness3*w_mult)/cos(dip_strike(3,1))); 
fthick4=abs((fault_thickness4*w_mult)/cos(dip_strike(4,1))); 
fthick5=abs((fault_thickness5*w_mult)/cos(dip_strike(5,1))); 
fthick6=abs((fault_thickness6*w_mult)/cos(dip_strike(6,1))); 
fthick7=abs((fault_thickness7*w_mult)/cos(dip_strike(7,1))); 
fthick8=abs((fault_thickness8*w_mult)/cos(dip_strike(8,1))); 
fthick9=abs((fault_thickness9*w_mult)/cos(dip_strike(9,1))); 
fthick10=abs((fault_thickness10*w_mult)/cos(dip_strike(10,1))); 
fthick11=abs((fault_thickness11*w_mult)/cos(dip_strike(11,1))); 

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
int3=-10000;

p11=[x01+cos(dip_strike(1,2)),y01+sin(dip_strike(1,2)),z01]; %point taken down strike
p12=[x02+cos(dip_strike(2,2)),y02+sin(dip_strike(2,2)),z02];
p13=[x03+cos(dip_strike(3,2)),y03+sin(dip_strike(3,2)),z03];
p14=[x04+cos(dip_strike(4,2)),y04+sin(dip_strike(4,2)),z04];
p15=[x05+cos(dip_strike(5,2)),y05+sin(dip_strike(5,2)),z05];
p16=[x06+cos(dip_strike(6,2)),y06+sin(dip_strike(6,2)),z06];
p17=[x07+cos(dip_strike(7,2)),y07+sin(dip_strike(7,2)),z07];
p18=[x08+cos(dip_strike(8,2)),y08+sin(dip_strike(8,2)),z08];
p19=[x09+cos(dip_strike(9,2)),y09+sin(dip_strike(9,2)),z09];
p110=[x010+cos(dip_strike(10,2)),y010+sin(dip_strike(10,2)),z010];
p111=[x011+cos(dip_strike(11,2)),y011+sin(dip_strike(11,2)),z011];

p21=[x01+cos(dip_strike(1,2)+(pi/2)),y01+sin(dip_strike(1,2)+(pi/2)),slope1+z01]; %point taken down dip
p22=[x02+cos(dip_strike(2,2)+(pi/2)),y02+sin(dip_strike(2,2)+(pi/2)),slope2+z02];
p23=[x03+cos(dip_strike(3,2)+(pi/2)),y03+sin(dip_strike(3,2)+(pi/2)),slope3+z03]; 
p24=[x04+cos(dip_strike(4,2)+(pi/2)),y04+sin(dip_strike(4,2)+(pi/2)),slope4+z04];
p25=[x05+cos(dip_strike(5,2)+(pi/2)),y05+sin(dip_strike(5,2)+(pi/2)),slope5+z05];
p26=[x06+cos(dip_strike(6,2)+(pi/2)),y06+sin(dip_strike(6,2)+(pi/2)),slope6+z06];
p27=[x07+cos(dip_strike(7,2)+(pi/2)),y07+sin(dip_strike(7,2)+(pi/2)),slope7+z07];
p28=[x08+cos(dip_strike(8,2)+(pi/2)),y08+sin(dip_strike(8,2)+(pi/2)),slope8+z08]; 
p29=[x09+cos(dip_strike(9,2)+(pi/2)),y09+sin(dip_strike(9,2)+(pi/2)),slope9+z09];
p210=[x010+cos(dip_strike(10,2)+(pi/2)),y010+sin(dip_strike(10,2)+(pi/2)),slope10+z010];
p211=[x011+cos(dip_strike(11,2)+(pi/2)),y011+sin(dip_strike(11,2)+(pi/2)),slope11+z011];

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
intnode=length(xc);
fault_planes(:,1)=((ortho1(1)*(xc-p01(1))+ortho1(2)*(yc-p01(2)))/ortho1(3))+z01;
zlt1b=fault_planes(:,1)+fthick1;
zlt1c=fault_planes(:,1)+fthick1/5;
zlt1d=fault_planes(:,1)+4*fthick1/5;
zlt1e=fault_planes(:,1)+2*fthick1/5;
zlt1f=fault_planes(:,1)+3*fthick1/5;

zlt2a=((ortho2(1)*(xc-p02(1))+ortho2(2)*(yc-p02(2)))/ortho2(3))+z02;
zlt2b=zlt2a+fthick2;
zlt2c=zlt1a+fthick2/5;
zlt2d=zlt1a+4*fthick2/5;
zlt2e=zlt1a+2*fthick2/5;
zlt2f=zlt1a+3*fthick2/5;

zlt3a=((ortho3(1)*(xc-p03(1))+ortho3(2)*(yc-p03(2)))/ortho3(3))+z03;
zlt3b=zlt3a+fthick3;
zlt3c=zlt1a+fthick3/5;
zlt3d=zlt1a+4*fthick3/5;
zlt3e=zlt1a+2*fthick3/5;
zlt3f=zlt1a+3*fthick3/5;

zlt4a=((ortho4(1)*(xc-p04(1))+ortho4(2)*(yc-p04(2)))/ortho4(3))+z04;
zlt4b=zlt4a+fthick4;
zlt4c=zlt1a+fthick4/5;
zlt4d=zlt1a+4*fthick4/5;
zlt4e=zlt1a+2*fthick4/5;
zlt4f=zlt1a+3*fthick4/5;

zlt5a=((ortho5(1)*(xc-p05(1))+ortho5(2)*(yc-p05(2)))/ortho5(3))+z05;
zlt5b=zlt5a+fthick5;
zlt5c=zlt1a+fthick5/5;
zlt5d=zlt1a+4*fthick5/5;
zlt5e=zlt1a+2*fthick5/5;
zlt5f=zlt1a+3*fthick5/5;

zlt6a=((ortho6(1)*(xc-p06(1))+ortho6(2)*(yc-p06(2)))/ortho6(3))+z06;
zlt6b=zlt6a+fthick6;
zlt6c=zlt1a+fthick6/5;
zlt6d=zlt1a+4*fthick6/5;
zlt6e=zlt1a+2*fthick6/5;
zlt6f=zlt1a+3*fthick6/5;

zlt7a=((ortho7(1)*(xc-p07(1))+ortho7(2)*(yc-p07(2)))/ortho7(3))+z07;
zlt7b=zlt7a+fthick7;
zlt7c=zlt1a+fthick7/5;
zlt7d=zlt1a+4*fthick7/5;
zlt7e=zlt1a+2*fthick7/5;
zlt7f=zlt1a+3*fthick7/5;

zlt8a=((ortho8(1)*(xc-p08(1))+ortho8(2)*(yc-p08(2)))/ortho8(3))+z08;
zlt8b=zlt8a+fthick8;
zlt8c=zlt1a+fthick8/5;
zlt8d=zlt1a+4*fthick8/5;
zlt8e=zlt1a+2*fthick8/5;
zlt8f=zlt1a+3*fthick8/5;

zlt9a=((ortho9(1)*(xc-p09(1))+ortho9(2)*(yc-p09(2)))/ortho9(3))+z09;
zlt9b=zlt9a+fthick9;
zlt9c=zlt1a+fthick9/5;
zlt9d=zlt1a+4*fthick9/5;
zlt9e=zlt1a+2*fthick9/5;
zlt9f=zlt1a+3*fthick9/5;

zlt10a=((ortho10(1)*(xc-p010(1))+ortho10(2)*(yc-p010(2)))/ortho10(3))+z010;
zlt10b=zlt10a+fthick10;
zlt10c=zlt1a+fthick10/5;
zlt10d=zlt1a+4*fthick10/5;
zlt10e=zlt1a+2*fthick10/5;
zlt10f=zlt1a+3*fthick10/5;

zlt11a=((ortho11(1)*(xc-p011(1))+ortho11(2)*(yc-p011(2)))/ortho11(3))+z011;
zlt11b=zlt11a+fthick11;
zlt11c=zlt1a+fthick11/5;
zlt11d=zlt1a+4*fthick11/5;
zlt11e=zlt1a+2*fthick11/5;
zlt11f=zlt1a+3*fthick11/5;

zltn=repmat(int3,length(xc),1);


thickness=zeros(maxnumlay, length(xc));
C=nan(maxnumlay,length(xc));
nlayers=zeros(1,length(xc));
depthdat=nan(maxnumlay+1,length(xc));
cdat=depthdat;

%omit any layer sections that would exist above the surface elevation
for i=1:length(xc);
    if zltn(i)>zci(i), zltn(i)=zci(i); end
    if zlt1a(i)>zci(i), zlt1a(i)=zci(i); end
    if zlt1b(i)>zci(i), zlt1b(i)=zci(i); end
    if zlt2a(i)>zci(i), zlt2a(i)=zci(i); end
    if zlt2b(i)>zci(i), zlt2b(i)=zci(i); end
    if zlt3a(i)>zci(i), zlt3a(i)=zci(i); end
    if zlt3b(i)>zci(i), zlt3b(i)=zci(i); end
    if zlt4a(i)>zci(i), zlt4a(i)=zci(i); end
    if zlt4b(i)>zci(i), zlt4b(i)=zci(i); end
    if zlt5a(i)>zci(i), zlt5a(i)=zci(i); end
    if zlt5b(i)>zci(i), zlt5b(i)=zci(i); end
    if zlt6a(i)>zci(i), zlt6a(i)=zci(i); end
    if zlt6b(i)>zci(i), zlt6b(i)=zci(i); end
    if zlt7a(i)>zci(i), zlt7a(i)=zci(i); end
    if zlt7b(i)>zci(i), zlt7b(i)=zci(i); end
    if zlt8a(i)>zci(i), zlt8a(i)=zci(i); end
    if zlt8b(i)>zci(i), zlt8b(i)=zci(i); end
    if zlt9a(i)>zci(i), zlt9a(i)=zci(i); end
    if zlt9b(i)>zci(i), zlt9b(i)=zci(i); end
    if zlt10a(i)>zci(i), zlt10a(i)=zci(i); end
    if zlt10b(i)>zci(i), zlt10b(i)=zci(i); end
    if zlt11a(i)>zci(i), zlt11a(i)=zci(i); end
    if zlt11b(i)>zci(i), zlt11b(i)=zci(i); end

    if zlt1c(i)>zci(i), zlt1c(i)=zci(i); end
    if zlt1d(i)>zci(i), zlt1d(i)=zci(i); end
    if zlt2c(i)>zci(i), zlt2c(i)=zci(i); end
    if zlt2d(i)>zci(i), zlt2d(i)=zci(i); end
    if zlt3c(i)>zci(i), zlt3c(i)=zci(i); end
    if zlt3d(i)>zci(i), zlt3d(i)=zci(i); end
    if zlt4c(i)>zci(i), zlt4c(i)=zci(i); end
    if zlt4d(i)>zci(i), zlt4d(i)=zci(i); end
    if zlt5c(i)>zci(i), zlt5c(i)=zci(i); end
    if zlt5d(i)>zci(i), zlt5d(i)=zci(i); end
    if zlt6c(i)>zci(i), zlt6c(i)=zci(i); end
    if zlt6d(i)>zci(i), zlt6d(i)=zci(i); end
    if zlt7c(i)>zci(i), zlt7c(i)=zci(i); end
    if zlt7d(i)>zci(i), zlt7d(i)=zci(i); end
    if zlt8c(i)>zci(i), zlt8c(i)=zci(i); end
    if zlt8d(i)>zci(i), zlt8d(i)=zci(i); end
    if zlt9c(i)>zci(i), zlt9c(i)=zci(i); end
    if zlt9d(i)>zci(i), zlt9d(i)=zci(i); end
    if zlt10c(i)>zci(i), zlt10c(i)=zci(i); end
    if zlt10d(i)>zci(i), zlt10d(i)=zci(i); end
    if zlt11c(i)>zci(i), zlt11c(i)=zci(i); end
    if zlt11d(i)>zci(i), zlt11d(i)=zci(i); end

    if zlt1e(i)>zci(i), zlt1e(i)=zci(i); end
    if zlt1f(i)>zci(i), zlt1f(i)=zci(i); end
    if zlt2e(i)>zci(i), zlt2e(i)=zci(i); end
    if zlt2f(i)>zci(i), zlt2f(i)=zci(i); end
    if zlt3e(i)>zci(i), zlt3e(i)=zci(i); end
    if zlt3f(i)>zci(i), zlt3f(i)=zci(i); end
    if zlt4e(i)>zci(i), zlt4e(i)=zci(i); end
    if zlt4f(i)>zci(i), zlt4f(i)=zci(i); end
    if zlt5e(i)>zci(i), zlt5e(i)=zci(i); end
    if zlt5f(i)>zci(i), zlt5f(i)=zci(i); end
    if zlt6e(i)>zci(i), zlt6e(i)=zci(i); end
    if zlt6f(i)>zci(i), zlt6f(i)=zci(i); end
    if zlt7e(i)>zci(i), zlt7e(i)=zci(i); end
    if zlt7f(i)>zci(i), zlt7f(i)=zci(i); end
    if zlt8e(i)>zci(i), zlt8e(i)=zci(i); end
    if zlt8f(i)>zci(i), zlt8f(i)=zci(i); end
    if zlt9e(i)>zci(i), zlt9e(i)=zci(i); end
    if zlt9f(i)>zci(i), zlt9f(i)=zci(i); end
    if zlt10e(i)>zci(i), zlt10e(i)=zci(i); end
    if zlt10f(i)>zci(i), zlt10f(i)=zci(i); end
    if zlt11e(i)>zci(i), zlt11e(i)=zci(i); end
    if zlt11f(i)>zci(i), zlt11f(i)=zci(i); end

    %     start sussing out layer data
    depthdat(1,i)=zci(i); cdat(1,i)=1;
    depthdat(2,i)=zlt1b(i); cdat(2,i)=0;
    depthdat(3,i)=zlt1d(i); cdat(3,i)=2;
    depthdat(4,i)=zlt1f(i); cdat(4,i)=4;
    depthdat(5,i)=zlt1e(i); cdat(5,i)=5;
    depthdat(6,i)=zlt1c(i); cdat(6,i)=3;    
    depthdat(7,i)=zlt1a(i); cdat(7,i)=1;
    if z02>-100000
        depthdat(8,i)=zlt2b(i); cdat(8,i)=0;
        depthdat(9,i)=zlt2d(i); cdat(9,i)=2;
        depthdat(10,i)=zlt2f(i); cdat(10,i)=4;
        depthdat(11,i)=zlt2e(i); cdat(11,i)=5;
        depthdat(12,i)=zlt2c(i); cdat(12,i)=3;
        depthdat(13,i)=zlt2a(i); cdat(13,i)=1;
    end
    if z03>-100000
        depthdat(14,i)=zlt3b(i); cdat(14,i)=0;
        depthdat(15,i)=zlt3d(i); cdat(15,i)=2;
        depthdat(16,i)=zlt3f(i); cdat(16,i)=4;
        depthdat(17,i)=zlt3e(i); cdat(17,i)=5;
        depthdat(18,i)=zlt3c(i); cdat(18,i)=3;
        depthdat(19,i)=zlt3a(i); cdat(19,i)=1;
    end
    if z04>-100000
        depthdat(20,i)=zlt4b(i); cdat(20,i)=0;
        depthdat(21,i)=zlt4d(i); cdat(21,i)=2;
        depthdat(22,i)=zlt4f(i); cdat(22,i)=4;
        depthdat(23,i)=zlt4e(i); cdat(23,i)=5;
        depthdat(24,i)=zlt4c(i); cdat(24,i)=3;
        depthdat(25,i)=zlt4a(i); cdat(25,i)=1;
    end
    if z05>-100000
        depthdat(26,i)=zlt5b(i); cdat(26,i)=0;
        depthdat(27,i)=zlt5d(i); cdat(27,i)=2;
        depthdat(28,i)=zlt5f(i); cdat(28,i)=4;
        depthdat(29,i)=zlt5e(i); cdat(29,i)=5;
        depthdat(30,i)=zlt5c(i); cdat(30,i)=3;
        depthdat(31,i)=zlt5a(i); cdat(31,i)=1;
    end
    if z06>-100000
        depthdat(32,i)=zlt6b(i); cdat(32,i)=0;
        depthdat(33,i)=zlt6d(i); cdat(33,i)=2;
        depthdat(34,i)=zlt6f(i); cdat(34,i)=4;
        depthdat(35,i)=zlt6e(i); cdat(35,i)=5;
        depthdat(36,i)=zlt6c(i); cdat(36,i)=3;
        depthdat(37,i)=zlt6a(i); cdat(37,i)=1;
    end
    if z07>-100000
        depthdat(38,i)=zlt7b(i); cdat(38,i)=0;
        depthdat(39,i)=zlt7d(i); cdat(39,i)=2;
        depthdat(40,i)=zlt7f(i); cdat(40,i)=4;
        depthdat(41,i)=zlt7e(i); cdat(41,i)=5;
        depthdat(42,i)=zlt7c(i); cdat(42,i)=3;
        depthdat(43,i)=zlt7a(i); cdat(43,i)=1;
    end
    if z08>-100000
        depthdat(44,i)=zlt8b(i); cdat(44,i)=0;
        depthdat(45,i)=zlt8d(i); cdat(45,i)=2;
        depthdat(46,i)=zlt8f(i); cdat(46,i)=4;
        depthdat(47,i)=zlt8e(i); cdat(47,i)=5;
        depthdat(48,i)=zlt8c(i); cdat(48,i)=3;
        depthdat(49,i)=zlt8a(i); cdat(49,i)=1;
    end
    if z09>-100000
        depthdat(50,i)=zlt9b(i); cdat(50,i)=0;
        depthdat(51,i)=zlt9d(i); cdat(51,i)=2;
        depthdat(52,i)=zlt9f(i); cdat(52,i)=4;
        depthdat(53,i)=zlt9e(i); cdat(53,i)=5;
        depthdat(54,i)=zlt9c(i); cdat(54,i)=3;
        depthdat(55,i)=zlt9a(i); cdat(55,i)=1;
    end
    if z010>-100000
        depthdat(56,i)=zlt10b(i); cdat(56,i)=0;
        depthdat(57,i)=zlt10d(i); cdat(57,i)=2;
        depthdat(58,i)=zlt10f(i); cdat(58,i)=4;
        depthdat(59,i)=zlt10e(i); cdat(59,i)=5;
        depthdat(60,i)=zlt10c(i); cdat(60,i)=3;
        depthdat(61,i)=zlt10a(i); cdat(61,i)=1;
    end
    if z011>-100000
        depthdat(62,i)=zlt11b(i); cdat(62,i)=0;
        depthdat(63,i)=zlt11d(i); cdat(63,i)=2;
        depthdat(64,i)=zlt11f(i); cdat(64,i)=4;
        depthdat(65,i)=zlt11e(i); cdat(65,i)=5;
        depthdat(66,i)=zlt11c(i); cdat(66,i)=3;
        depthdat(67,i)=zlt11a(i); cdat(67,i)=1;
    end
    
    depthdat(68,i)=zltn(i); cdat(68,i)=0;
    
    [Y,I]=sort(depthdat(:,i),'descend');
    depthdat(:,i)=depthdat(I,i);
    cdat(:,i)=cdat(I,i);
    ddat= depthdat(~isnan(depthdat(:,i)),i);
    ccdat= cdat(~isnan(cdat(:,i)),i);
    n=2; tic=1; tally=0; sally=tally; bally=sally;
    while n <= length(ddat)
        if ccdat(n)==0 %top of outer fault layer, bottom of host layer
            if tally==0 
                thickness(tic,i)=ddat(n-1)-ddat(n);
                C(tic,i)=c1;
                tic=tic+1;
            end
            if tally>0 && sally==0
                thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
                C(tic,i)=c2;
            end
            tally=tally+1;
            n=n+1; continue
        end
        
        if ccdat(n)==2 %top of middle fault layer, bottom of outer fault layer
            if sally==0 %terminate the outer layer thickness
            thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
            C(tic,i)=c2;
            tic=tic+1;
            end
            if sally>0 && bally==0 % accumulate the middle layer thicknesses
                thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
                C(tic,i)=c3;
            end
            sally=sally+1;
            tally=0;
            n=n+1; continue
        end
        
        if ccdat(n)==4 %top of core fault layer, bottom of middle fault layer
            if bally==0 %terminate the middle layer thickness
            thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
            C(tic,i)=c3;
            tic=tic+1;
            end
            if bally>0 % accumulate the core layer thicknesses
                thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
                C(tic,i)=c4;
            end
            bally=bally+1;
            tally=0; sally=tally;
            n=n+1; continue
        end
%         layer bottoms
        if ccdat(n)==5 %bottom of core layer, top of middle layer
            bally=bally-1;
            thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
            C(tic,i)=c4;
            n=n+1;
            if bally==0
                tic=tic+1;
            end
        end
        if ccdat(n)==3 %bottom of middle layer, top of outer layer
            sally=sally-1;
            thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
            C(tic,i)=c3;
            n=n+1;
            if sally==0
                tic=tic+1;
            end
        end
        
        if ccdat(n)==1 %bottom of outer layer, top of host
            tally=tally-1;
            thickness(tic,i)=thickness(tic,i)+ddat(n-1)-ddat(n);
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

% clip layers for plot
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


