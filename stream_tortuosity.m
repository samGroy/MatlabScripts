function c=stream_tortuosity(basenm,ts,wavelength,x,y,overprint)

% SGR 6/23/2013
% Plot data for a stream profile. The profile follows the downstream path
% from a chosen point in the domain (x,y coords). Returns stream profile,
% tortuosity with respect to a chosen wavelength, discharge, and channel
% width calculated both by Finnegan (2005) and the empirical relation to
% flux. The channel path is also highlighted in 3D space.

% Find and open data files
filesys=[''];
% ts=5;
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.nodes' ];
nodfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
filenm= [filesys basenm '.tri' ];
tfid=fopen(filenm,'r');
if tfid<=0, fclose(nfid); error('Unable to open triangle file'),end
fprintf('CPLOTAREA: Reading data ...\n');

for i=1:ts
  tm=fscanf(netfid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm);
  tm=fscanf(nodfid,'%f',1);
  tm=fscanf(zfid,'%f',1);
  tm = fscanf(tfid,'%f',1);
  nt = fscanf(tfid,'%d',1);
  t=fscanf(tfid,'%f',[9,nt]); 
  intnodes= fscanf(netfid,'%d',1);
  allnodes= fscanf(zfid,'%d',1);
  allnodes= fscanf(nodfid,'%d',1); 
  netdat=fscanf(netfid,'%d',[1,intnodes]);
  nodedat=fscanf(nodfid,'%f',[4,allnodes]);
  zdat=fscanf(zfid,'%f',[1,allnodes]);
end
% If it isn't used, turn off the option to overprint any proceeding
% tortuosity calculations
if nargin == 5, overprint = 0; end
% Convert the triangle data used to plot the base topo map
tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
zdat=zdat';
zorf=zeros(allnodes,1);

[id,steparc, nodedat, netdat]=tortuosityprofile5(basenm,ts,wavelength,x,y);
% [id2,steparc2, nodedat2, netdat2]=tortuosityprofile4(basenm,ts,wavelength2,x,y);
if overprint == 0
    fhandle2=figure;%('Position',[100,600,1000,1000]);
    h=trisurf(tri,nodedat(:,1),nodedat(:,2),zorf,zdat'); hold on;
    title(['stream path, ',num2str(wavelength),' m, ',num2str(tm),' years']);
    axis([min(nodedat(:,1)) max(nodedat(:,1)) min(nodedat(:,2)) max(nodedat(:,2)) min(zdat) (max(zdat))]);
    shading interp;
    view([-44 50]);
    daspect([1 1 100]);
end
% h=plot3(nodedat(id(:,1),1),nodedat(id(:,1),2),zdat(id(:,1)),'b','LineWidth',10); 
%k=scatter3(nodedat(id(:,1),1),nodedat(id(:,1),2),zdat(id(:,1)),(empwidth(id(:,1))+10),steparc(:,3),'filled'); 
if overprint ~= 0, figure(overprint); end
% k=scatter3(nodedat(id(:,1),1),nodedat(id(:,1),2),zdat(id(:,1)),20,steparc(:,3),'filled'); hold on;
k=scatter3(steparc(:,4),steparc(:,5),zdat(steparc(:,6)),7,steparc(:,3),'filled'); hold on; %use 7 for regimetests, 5 for EHS, apparently rounds to 5's
colormap([gray(64);jet(64)])
% Initially, both CDatas are equal to Z.
m = 64; % 64-elements is each colormap
% cmin = min(steparc(:,3)); %Use these two if you want the script to choose color scale limits based on the min and max of data automatically.
% cmax = max(steparc(:,3));
cmin = 0; %Use these two if you want to set the color scale limites for tortuosity, useful fr comparison
cmax = 0.5;
if overprint == 0 % Only define these limits if you're making the surface map
    zmin = min(zdat); zmax=max(zdat);
    % CData for surface
    C1 = min(m,round((m-1)*(zdat-zmin)/(zmax-zmin))+1);
end
% CData for tortuosity
C2 = 64+min(m,round((m-1)*(steparc(:,3)-cmin)/(cmax-cmin))+1);


% Update the CDatas for each object.
if overprint == 0
    set(h,'CData',C1);
end
set(k,'CData',C2);

% Change the CLim property of axes so that it spans the 
% CDatas of both objects.
% caxis([min(C1(:)) max(C2(:))])
if overprint == 0
    caxis([min(C1(:)) 128])
end

mean(steparc(:,3))
2*std(steparc(:,3))
fclose(netfid);
fclose(nodfid);
fclose(zfid);
