function [id,tort,trig,nodedat,netdat]=tortuosityprofile(basenm,ts,wavelength,x,y)
% basenm='fractureset01';
% x=5000;y=5000;
% wavelength=100;
filesys=[''];
% ts=5;

% collect stream network data
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.nodes' ];
nodfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
fprintf('CPLOTAREA: Reading data ...\n');

for i=1:ts
  tm=fscanf(netfid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm);
  tm=fscanf(nodfid,'%f',1);
  tm=fscanf(zfid,'%f',1);
  intnodes= fscanf(netfid,'%d',1);
  allnodes= fscanf(zfid,'%d',1);
  allnodes= fscanf(nodfid,'%d',1); 
  netdat=fscanf(netfid,'%d',[1,intnodes]);
  nodedat=fscanf(nodfid,'%f',[4,allnodes]);
  zdat=fscanf(zfid,'%f',[1,allnodes]);
end

nodedat=nodedat';
netdat=(netdat+1)';

% find the closest interior node to the input xy coordinates to start the profiles
for iiii=1:allnodes % collect interior node data
    if nodedat(iiii,4)==0, intnodedat(iiii,:)=nodedat(iiii,:); end
end
dx=intnodedat(:,1)-x; dx=dx.^2;
dy=intnodedat(:,2)-y; dy=dy.^2;
dist1=sqrt(dx+dy);
[a,b]=min(dist1); 
secperyear=31556926;
point_in_interior=1;
o=1;
%sets up a vector with nearest downslope IDs
while point_in_interior
    id(o,1)=b; % the id vector keeps track of downstream node order
    if o==1, id(o,2)=0; id(o,3)=0; end
    o=o+1;
    if nodedat(b,4)~=0,point_in_interior=0;break % allows inclusion of
%     boundary outlet node in id list if located here
    end
    dx=nodedat(netdat(b),1)-nodedat(b,1); dx=dx.^2;
    dy=nodedat(netdat(b),2)-nodedat(b,2); dy=dy.^2;
    id(o,2)=sqrt(dx+dy);
    id(o,3)=id(o,2)+id(o-1,3);
b=netdat(b);
%     if nodedat(b,4)~=0,point_in_interior=0;break 
end

trig=nan(intnodes,2);
tort=nan(length(id),3);
% finds angle between node and adjacent downstream
% node, angle increases ccw to 360/0 due east
for ii=1:intnodes
   if nodedat(ii,1)<nodedat(netdat(ii),1) && nodedat(ii,2)<nodedat(netdat(ii),2) % quadrant I
       trig(ii)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))));
   end
   if nodedat(ii,1)>nodedat(netdat(ii),1) && nodedat(ii,2)<nodedat(netdat(ii),2) % quadrant II
       trig(ii)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))))+pi/2;
   end
   if nodedat(ii,1)>nodedat(netdat(ii),1) && nodedat(ii,2)>nodedat(netdat(ii),2) % quadrant III
       trig(ii)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))))+2*pi/2;
   end
   if nodedat(ii,1)<nodedat(netdat(ii),1) && nodedat(ii,2)>nodedat(netdat(ii),2) % quadrant IV
       trig(ii)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))))+3*pi/2;
   end
      trig(ii,2)=((nodedat(netdat(ii),2)-nodedat(ii,2)).^2+(nodedat(netdat(ii),1)-nodedat(ii,1)).^2).^.5; %measures node to node distance

%    if trig(ii,3)>intnodes
%        trig(ii,4)=trig(ii,3);
%    else
%    trig(ii,4)=netdat(netdat(ii));
%    end
% try to replace above with a vector collecting downstream IDs and
% euclidean distances
end

% find change in orientation from upstream node pairs and downstream node
% pairs and euclidean distance between the first upstream node and the
% first downstream node. Breaks when the downstream outlet node is reached.
for i=1:length(id(:,1))
    if i==1
        if mean(trig(:,2))-mean(trig(:,2))/4>wavelength
            wavelength=mean(trig(:,2)); 
            fprintf('Chosen wavelength is too far below average node spacing.\nSwitching wavelength to average spacing distance: %f \n',mean(trig(:,2)));
        end
    end
    dx=nodedat(id(i,1),1)-nodedat(id(:,1),1); dx=dx.^2;
    dy=nodedat(id(i,1),2)-nodedat(id(:,1),2); dy=dy.^2;
    dist2=zeros(length(id(:,1)),2);
    dist2(:,1)=sqrt(dx+dy);
    dist2(i:end,2)=abs(wavelength-dist2(i:end,1));
    [e,f]=min(dist2(i+1:end,2)); f=f+i;
    if dist2(f)<wavelength/3; break % if distance to next point is less than 3/4 the wavelength distance, break.
    end
    if nodedat(id(f,1),4)~=0, break
    end
tort(i,1)=abs(trig(id(i,1))-trig(id(f,1)));
    tort(i,2)=dist2(f,1);
    tort(i,3)=tort(i,1)/tort(i,2);
%     if i==1 %this makes tort(:,2) equal the cumulate distance from the starting point, this is used for the plot
%         tort(i,2)=tort(i,2); % zero total distance for first node
%     else
%         tort(i,2)=tort(i,2)+tort(i-1,2); % total distance used for tortuosity plot
%     end
end



fclose(netfid);
fclose(nodfid);
fclose(zfid);
%fclose(afid);
%if wfid>0, fclose(wfid); end
