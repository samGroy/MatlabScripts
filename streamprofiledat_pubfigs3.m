function c=streamprofiledat_pubfigs2(basenm,ts,wavelength1,wavelength2,x,y)

% SGR 8/2012
% Plot data for a stream profile. The profile follows the downstream path
% from a chosen point in the domain (x,y coords). Returns stream profile,
% tortuosity with respect to a chosen wavelength, discharge, and channel
% width calculated both by Finnegan (2005) and the empirical relation to
% flux. The channel path is also highlighted in 3D space.

% basenm='fractureset01';
% x=10000;y=10000;
% wavelength=500;
filesys=[''];
% ts=5;
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.nodes' ];
nodfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
filenm= [filesys basenm '.area' ];
afid=fopen(filenm,'r');
filenm= [filesys basenm '.chanwid' ];
wfid=fopen(filenm,'r');
filenm= [filesys basenm '.q' ];
qfid=fopen(filenm,'r');
filenm= [filesys basenm '.tri' ];
tfid=fopen(filenm,'r');
if tfid<=0,  error('Unable to open triangle file'),end
fprintf('CPLOTAREA: Reading data ...\n');

for i=1:ts
  tm=fscanf(netfid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm);
  tm=fscanf(nodfid,'%f',1);
  tm=fscanf(zfid,'%f',1);
  tm=fscanf(afid,'%f',1);
  tm=fscanf(qfid,'%f',1);
  tm = fscanf(tfid,'%f',1);
  nt = fscanf(tfid,'%d',1);
  t=fscanf(tfid,'%f',[9,nt]); 
  intnodes= fscanf(afid,'%d',1);
  allnodes= fscanf(qfid,'%d',1);
  intnodes= fscanf(netfid,'%d',1);
  allnodes= fscanf(zfid,'%d',1);
  allnodes= fscanf(nodfid,'%d',1); 
  netdat=fscanf(netfid,'%d',[1,intnodes]);
  nodedat=fscanf(nodfid,'%f',[4,allnodes]);
  zdat=fscanf(zfid,'%f',[1,allnodes]);
  area=fscanf(afid,'%f',[1,intnodes]);
  if wfid<0
      width=nan(1,allnodes); end
  if wfid>0
      tm=fscanf(wfid,'%f',1);
      allnodes= fscanf(wfid,'%d',1);
      width=fscanf(wfid,'%f',[1,allnodes]);
  end
  q=fscanf(qfid,'%f',[1,allnodes]);
end
secperyear=31556926;
empwidth=10*(q/secperyear).^.5;
tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
zdat=zdat';
zorf=zeros(allnodes,1);

[id,steparc, nodedat, netdat]=tortuosityprofile4(basenm,ts,wavelength1,x,y);
[id2,steparc2, nodedat2, netdat2]=tortuosityprofile4(basenm,ts,wavelength2,x,y);
fhandle2=figure;%('Position',[100,600,1000,1000]);
h=trisurf(tri,nodedat(:,1),nodedat(:,2),zorf,zdat'); hold on; 
title(['stream path, ',num2str(tm),' years']);
axis([min(nodedat(:,1)) max(nodedat(:,1)) min(nodedat(:,2)) max(nodedat(:,2)) min(zdat)-1000 (max(zdat)+1000)]);
shading flat;
view([0 90]);
% h=plot3(nodedat(id(:,1),1),nodedat(id(:,1),2),zdat(id(:,1)),'b','LineWidth',10); 
% k=scatter3(nodedat(id(:,1),1),nodedat(id(:,1),2),zdat(id(:,1)),(empwidth(id(:,1))+10),steparc(:,3),'filled'); 
% k=scatter3(nodedat(id(:,1),1),nodedat(id(:,1),2),5+zorf(id(:,1)),(empwidth(id(:,1))+10),steparc(:,3),'filled'); 
k=plot3(nodedat(id(:,1),1),nodedat(id(:,1),2),5+zorf(id(:,1)),'w','LineWidth',2); 
%j=plot3(idticker(:,1),idticker(:,2),10+zorf(1:length(idticker)),'o','MarkerEdgeColor','k','MarkerFaceColor','w','MarkerSize',4);
colormap([gray(64);hot(64)])
m = 128; % 64-elements is each colormap

cmin = min(steparc(:,3));
cmax = max(steparc(:,3));
zmin = min(zdat); zmax=max(zdat);
% CData for surface
C1 = min(m,round((m-1)*(zdat-zmin)/(zmax-zmin))+1); 
% CData for pcolor
C2 = 64+min(m,round((m-1)*(steparc(:,3)-cmin)/(cmax-cmin))+1);


% Update the CDatas for each object.
set(h,'CData',C1);
set(k,'CData',C2);


% Change the CLim property of axes so that it spans the 
% CDatas of both objects.
caxis([min(C1(:)) max(C2(:))])
hold off;

%%%%%%%%%%%%%%%%%%%%%%%% Currently disabled for editing purposes SGR 2/27/13
% Initially, both CDatas are equal to Z.
% fhandle2b=figure;%('Position',[100,600,1000,1000]);
% hb=trisurf(tri,nodedat2(:,1),nodedat2(:,2),zorf,zdat'); hold on; 
% title(['stream path, ',num2str(tm),' years']);
% axis([min(nodedat2(:,1)) max(nodedat2(:,1)) min(nodedat2(:,2)) max(nodedat2(:,2)) min(zdat)-1000 (max(zdat)+1000)]);
% shading flat;
% view([0 90]);
% % h=plot3(nodedat(id(:,1),1),nodedat(id(:,1),2),zdat(id(:,1)),'b','LineWidth',10); 
% kb=scatter3(nodedat2(id2(:,1),1),nodedat2(id2(:,1),2),zdat(id2(:,1)),(empwidth(id2(:,1))+10),steparc2(:,3),'filled'); 
% colormap([gray(64);jet(64)])
% m = 64; % 64-elements is each colormap
% 
% cmin = min(steparc2(:,3));
% cmax = max(steparc2(:,3));
% zmin = min(zdat); zmax=max(zdat);
% % CData for surface
% C1 = min(m,round((m-1)*(zdat-zmin)/(zmax-zmin))+1); 
% % CData for pcolor
% C2 = 64+min(m,round((m-1)*(steparc2(:,3)-cmin)/(cmax-cmin))+1);
% 
% 
% % Update the CDatas for each object.
% set(hb,'CData',C1);
% set(kb,'CData',C2);
% 
% 
% % Change the CLim property of axes so that it spans the 
% % CDatas of both objects.
% caxis([min(C1(:)) max(C2(:))])
% hold off;

fhandle3=figure;%('Position', [100,100,350,350]);
subplot(3,1,1), [ax,h1a,h2a]=plotyy(id(:,3),steparc(:,3),id(:,3),zdat(id(:,1)),'plot'); hold on;

subplot(3,1,1), [ax,h1b,h2b]=plotyy(id2(:,3),steparc2(:,3),id(:,3),zdat(id(:,1)),'plot'); hold on;
title(['tortuosity, ',num2str(wavelength1),' m (dotted) and',num2str(wavelength2),' m (solid) wavelengths']); 
set(h1a,'Color','k','LineStyle',':'); 
set(h1b,'Color','k','LineStyle','-'); 
set(get(ax(1),'Ylabel'),'String','Tortuosity') 



hold off;

subplot(3,1,2), [ax,h3,h4]=plotyy(id(1:end-1,3),area(id(1:end-1,1)),id(1:end-1,3),zdat(id(1:end-1,1)),'plot'); hold on;
% title('drainage area');
set(h3,'Color','k');
set(h4,'Color','b'); hold off;
set(get(ax(1),'Ylabel'),'String','Drainage Area (m.^2)') 
set(get(ax(2),'Ylabel'),'String','Elevation') 

% if wfid>0
% subplot(3,1,3), [ax,h5,h6]=plotyy(id(:,3),width(id(:,1)),id(:,3),zdat(id(:,1)),'plot'); hold on;
% set(h5,'Color','g');
% set(h6,'Color','b'); 
% end
% 

subplot(3,1,3), [ax,h7,h8]=plotyy(id(1:end-1,3),empwidth(id(1:end-1,1)),id(1:end-1,3),zdat(id(1:end-1,1)),'plot'); hold on;
% title('channel width, Finnegan(g) and Empirical(r)');
set(h7,'Color','k');
set(h8,'Color','b'); hold off;
set(get(ax(1),'Ylabel'),'String','Channel Width (m)') 



fclose(netfid);
fclose(nodfid);
fclose(zfid);
fclose(afid);
if wfid>0, fclose(wfid); end
