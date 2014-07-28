function [knickrate,knick,flux]=xsect2(basenm,total_ts,dip)
% SGR 2/2013, function used to 1) create profiles along strike, 2)
% calculate knickpoint migration rates, and 3) calculate flux at the outlet
% node. Used to compare channels with varying fault dip.
% find average node spacing, x limits, y limits, total time steps

% SGR 8/2013 DOESN'T WORK if lakefill = 0 and tlim is used. the max flux
% and upstream algo get cut off where flow is cut off by sediment buildup.
data = textread(['' basenm '.inputs'], '%s', 'delimiter', '\n');
step=str2double(data(27,1));
knick=[];cross=[];srt=[];flux=[];
xlim=str2double(data(21,1));
ylim=str2double(data(23,1));
 if nargin==1, total_ts=str2double(data(5,1))/str2double(data(7,1)); end
knick(1:total_ts+1,3)=[0:total_ts]*str2double(data(7,1));
nfile=fopen(['' basenm '.nodes'],'r');
zfile=fopen(['' basenm '.z'],'r');
qfile=fopen(['' basenm '.q'],'r');
netfile=fopen(['' basenm '.net'],'r');

thick=5;
depth_to_footwall=100-((thick*5)/cos(deg2rad(dip)))*3/5;
figure(10); if depth_to_footwall>0, plot([0 ylim],[depth_to_footwall depth_to_footwall],':k'); hold on; end
set(figure(10),'Position',[500 500 1000 300],'Color','white')
xlabel('Distance from outlet (m)','FontSize',12);
ylabel('elevation (m)','FontSize',12);
% chosen x coordinate
% x_in=970;%1025(30dip);1695(noF);1015(>30dip);1020;1045(20dip);1080(10dip);1170(5dip);1055(15dip)
x_in=round((xlim/2)+(thick*3)/sin(deg2rad(dip))); %find general position of outlet, not actually used to search for fault controlled channel
%x_in=426;
y_in=ylim-step;
if x_in>xlim
    fprintf('x input is greater than x limit %f\n Try lower, you fool.\n',xlim)
end
tic=1;
for time=0:1:total_ts
    cross=[];
    fprintf('time step %f...\n',time);
    t=fscanf(nfile,'%f',1); t=fscanf(zfile,'%f',1); t=fscanf(qfile,'%f',1); t=fscanf(netfile,'%f',1);
    nodecount=fscanf(nfile,'%f',1); nodecount=fscanf(zfile,'%f',1); nodecount=fscanf(qfile,'%f',1); inodecount=fscanf(netfile,'%f',1);
    n=(fscanf(nfile,'%f',[4,nodecount]))';
    z=(fscanf(zfile,'%f',[1,nodecount]))';
    q=(fscanf(qfile,'%f',[1,nodecount]))';
    net=(fscanf(netfile,'%f',[1,inodecount]))';
    
    %     if time==0
    %         dx=n(:,1)-x_in; dx=dx.^2;
    %         dy=n(:,2); dy=dy.^2;
    %         dist=sqrt(dx+dy);
    %         [a,b]=min(dist); %index of closest starting node
    %         i=1;
    %     else
    [a,b]=max(q); %index of greatest discharge rate, indicates fault controlled outlet
%     if time == 1
%         bold = b;
%     elseif time > 1
%         b = bold;
%     end
    i=2;
    %     end
    %     while tic
    %         if z(net(b))<100,break
    %         else
    %             ba(i,1:4)=[b z(b) n(b,1) n(b,2)]; i=i+1;
    %             b=net(b);
%         end
%     end
%     find start node using least squares
    while tic
        if i==2; %indicate boundary node, which has no flux
            cross(i-1,1:2)=n(net(b),1:2);
            cross(i-1,3)=z(net(b));
            cross(i-1,4)=q(net(b));
        end
        cross(i,1:2)=n(b,1:2);
        cross(i,3)=z(b);
        cross(i,4)=q(b);
        if z(b)==100, break, end
        b=find(net==b-1);
%         [a,yy]=max(n(b,2)); 
        [a,qq]=max(q(b)); 
%         if n(b,2)<cross(i,2)
        b=b(qq);
%         else
%             b=b(yy);
%         end
        if isempty(b) || (n(b,2)<cross(i,2) && cross(i,2)>990), break, end
        i=i+1;
    end
 
    cross(2:end,5)=diff(cross(:,3))./diff(cross(:,2)); %differential for slope
    cross(2:end,6)=diff(cross(:,5))./(diff(cross(:,2))); %differential for curvature
    %     cross(2:end,6)=(diff(cross(:,5))./(diff(cross(:,2))));
    %     [val,index2]=min(cross(:,5));
%     crossmooth=cross(:,4);
%     for aa=11:1:length(cross(:,3))-11
%         crossmooth(aa,1)=mean(cross(aa-10:aa+10,3));
%     end
%     srt=sort(crossmooth(:,1),'descend');
%     index2=find(crossmooth(:,1)==srt(1));
    [b,index2]=min(abs(depth_to_footwall-cross(:,3)));
    knick(tic,1)=cross(index2,2);
%     index3=find(cross(:,5)>1E-4,1,'last');
    index3=find(cross(:,3)==100,1,'first');
    if isempty(index3)
        knick(tic,2)=0;
    else
        knick(tic,2)=cross(index3,2);
    end
    flux(tic,1)=cross(2,4);
%     if time == 5 || time==10 || time==15 || time == 20 || time==25 || time == 30 || time == 40 || time == 50
%      if time==10 || time == 20 || time == 30 || time == 40 || time == 50
        figure(10)
%         [ax,h1a,h2a]=plotyy(cross(1:index3,2),cross(1:index3,3),cross(1:index3,2),cross(1:index3,4));hold on;
%         set(h1a,'Color','b','LineWidth',2);
%         set(h2a,'Color','b','LineWidth',2,'LineStyle',':');
        plot(cross(:,2),cross(:,3),'k','LineWidth',1); hold on;
        plot(knick(tic,1),cross(index2,3),'or')
%         plot(cross(:,2),cross(:,4),'r','LineWidth',2); hold on;
%         plot(cross(:,2),cross(:,5),'k','LineWidth',2); hold on;
        %         plot(cross(:,2),cross(:,6),'k','LineWidth',2); hold on;
        
%     end
    tic=tic+1;
end
index4=find(knick(:,2)<990 & knick(:,2)>0,1,'last');
index5=find(diff(knick(2:end,1))>50 | diff(knick(2:end,1))<0,1,'first');
if isempty(index5),index5=index4; end
if index5==index4-1, index5=index4; end
knickrate(1:2,1)=[mean(diff(knick(2:index5,1))./diff(knick(2:index5,3))) mean(diff(knick(2:index4,2))./diff(knick(2:index4,3)))];
hold off;