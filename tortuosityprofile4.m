function [id,steparc,nodedat,netdat]=tortuosityprofile4(basenm,ts,wavelength,x,y)
% SGR update 10/31/2012 (Happy Halloween)
% changed algo to measure ratio of step length/arc length rather than track
% orientation change. Might improve analysis when switching wavelengths.
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
for iiii=1:allnodes % collect interior node data, 6/13 MIGHT BE ABLE TO SPEED UP WITH LOGIC
    if nodedat(iiii,4)==0, intnodedat(iiii,:)=nodedat(iiii,:); end
end
dx=intnodedat(:,1)-x; dx=dx.^2;
dy=intnodedat(:,2)-y; dy=dy.^2;
dist1=sqrt(dx+dy);
[a,b]=min(dist1); %index of closest starting node
secperyear=31556926;
point_in_interior=1;
o=1; ticker=1000;d=1; % ticker makes a mark on the tortuosity path so you can measure distance. Those data are put in id2. Default is 100.

%Define the stream path used to measure tortuosity:
%sets up a vector with nearest downslope IDs
while point_in_interior
    id(o,1)=b; % the id vector keeps track of downstream node order
    if o==1, id(o,2)=0; id(o,3)=0; end
    id(o,4)=nodedat(b,1); id(o,5)=nodedat(b,2); %keeps track of x, y coordinates of this node
    if id(o,3) < ticker+1.5 && id(o,3) > ticker-1.5 %if the stream path has traveled a specific distance (100 m intervals)
        id2(d,1:2)=id(o,4:5); %make a new matrix keeping track of the xy coordinates so you can plot tick marks on the stream path
        ticker=ticker+100; d=ticker/100;
    end
    o=o+1;
    if nodedat(b,4)~=0,point_in_interior=0;break % allows inclusion of
        %     boundary outlet node in id list if located here
    end
    dx=nodedat(netdat(b),1)-nodedat(b,1); dx=dx.^2; %least squares method of finding closest node
    dy=nodedat(netdat(b),2)-nodedat(b,2); dy=dy.^2;
    id(o,2)=sqrt(dx+dy); % find euclidean distance between this node and next node
    id(o,3)=id(o,2)+id(o-1,3); % find total distance to this node
    b=netdat(b); %the next downstream node index
    %     if nodedat(b,4)~=0,point_in_interior=0;break
end

%Calculate tortuosity for each interior point on the stream path
% find change in orientation from upstream node pairs and downstream node
% pairs and euclidean distance between the first upstream node and the
% first downstream node. Breaks when the downstream outlet node is reached.
steparc=zeros(length(id),5);
back=0; g=1; front=0;
[e,g]=min(abs(id(:,3)-wavelength/2)); %find first node 1/2 wavelength from starting node
% for i=1:f-i:length(id(:,1))
figure
while front<=length(id(:,1))
    if back==0
        if mean(id(2:end,2))/4>wavelength
            wavelength=mean(id(2:end,2)); 
            fprintf('Chosen wavelength is too far below average node spacing, you fool.\nSwitching wavelength to average spacing distance: %f \n',mean(id(2:end,2)));
        end
    end
    dx=nodedat(id(g,1),1)-nodedat(id(:,1),1); dx=dx.^2;
    dy=nodedat(id(g,1),2)-nodedat(id(:,1),2); dy=dy.^2;
    dist2=zeros(length(id(:,1)),2);
    dist2(:,1)=sqrt(dx+dy);
    dist2(:,2)=abs(wavelength/2-dist2(:,1)); %get distances from current node to all other nodes in stream network
    [e,back]=min(dist2(1:g-1,2));
    dx=nodedat(id(back,1),1)-nodedat(id(:,1),1); dx=dx.^2;
    dy=nodedat(id(back,1),2)-nodedat(id(:,1),2); dy=dy.^2;
    dist2=zeros(length(id(:,1)),2);
    dist2(:,1)=sqrt(dx+dy);
    dist2(:,2)=abs(wavelength-dist2(:,1)); %get distances from current node to all other nodes in stream network
    [e,front]=min(dist2(back:end,2));
    front=front+back; %find next node that hs the best fit distance from the current node
    if front>=length(id(:,1)), break, end
    dx2=nodedat(id(back,1),1)-nodedat(id(front,1),1); dx2=dx2.^2;
    dy2=nodedat(id(back,1),2)-nodedat(id(front,1),2); dy2=dy2.^2;
    steparc(g,1)=sqrt(dx2+dy2);
    steparc(g,2)=id(front,3)-id(back,3);
    steparc(g,3)=1-steparc(g,1)/steparc(g,2);
    steparc(g,4)=nodedat(id(g,1),1); steparc(g,5)=nodedat(id(i,1),2);
    xlim=[min(id(1:front,4)) max(id(1:front,4))]; ylim=[min(id(1:front,5)) max(id(1:front,5))];
    plot(id(back:front,4),id(back:front,5),'-b',id(g,4),id(g,5),'ok'); hold on;
    plot([nodedat(id(back,1),1) nodedat(id(front,1),1)],[nodedat(id(back,1),2) nodedat(id(front,1),2)],'-r');
    
%     if dist2(f)<wavelength/3; break % if distance to next point is less than 3/4 the wavelength distance, break.
%     end
    if nodedat(id(front,1),4)~=0, steparc(front,:)=steparc(front-1,:);break
    end
    g=g+1;
%     if i==1 %this makes tort(:,2) equal the cumulate distance from the starting point, this is used for the plot
%         tort(i,2)=tort(i,2); % zero total distance for first node
%     else
%         tort(i,2)=tort(i,2)+tort(i-1,2); % total distance used for tortuosity plot
%     end
end
title(['tortuosity iterations, ',num2str(wavelength),' m wavelength']);

fclose(netfid);
fclose(nodfid);
fclose(zfid);
%fclose(afid);
%if wfid>0, fclose(wfid); end
