function a=waveorient2(basenm,ts)
filesys=[''];
%ts=5;
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.nodes' ];
nodfid=fopen(filenm,'r');
filenm= [filesys basenm '.q' ];
qfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
% filenm= [filesys basenm '.tri' ];
% tfid=fopen(filenm,'r');

%nd(:,3)=nd(:,3)+5;
fprintf('CPLOTAREA: Reading data ...\n');

for i=1:ts
  tm=fscanf(netfid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm);
  tm=fscanf(nodfid,'%f',1);
  intnodes= fscanf(netfid,'%d',1);
  allnodes= fscanf(nodfid,'%d',1); 
  netdat=fscanf(netfid,'%d',[1,intnodes]);
  nodedat=fscanf(nodfid,'%f',[4,allnodes]);
  tm=fscanf(qfid,'%f',1);
  allnodes= fscanf(qfid,'%d',1); 
  q=fscanf(qfid,'%f',[1,allnodes]);

  
  tm=fscanf(zfid,'%f',1);
  allnodes= fscanf(zfid,'%d',1);
  z=fscanf(zfid,'%f',[1,allnodes]);
%   nt = fscanf(tfid,'%d',1);
%   t=fscanf(tfid,'%f',[9,nt]); 

end

% tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
nodedat=nodedat';
netdat=(netdat+1)';
z=z';
tally=1;
% for i=2:4:length(q)
%     if q(i)>=500
%         dat(tally,1)=nodedat(i,1);
%         dat(tally,2)=nodedat(i,2);
%         dat(tally,3)=i;
%         dat(tally,4)=q(i);
%         tally=tally+1;
%     end
% end
tally=1; tic=tally; toc=tic;
dx=diff(nodedat(:,1)).^2;
dy=diff(nodedat(:,2)).^2;
dist=sqrt(dx+dy);
minwave=min(dist);
maxwave=1000;
thetamax=179;
w(1:maxwave/10,1:thetamax)=0;
count(1:maxwave,1:thetamax)=0;
for i=1:length(nodedat)
    for theta=0:2:thetamax
        for wave=minwave:10:maxwave
            x=cos(deg2rad(theta+1))*wave+nodedat(i,1); if x<0, continue, end
            y=sin(deg2rad(theta+1))*wave+nodedat(i,2); if y<0, continue, end
            dx=nodedat(:,1)-x; dx=dx.^2;
            dy=nodedat(:,2)-y; dy=dy.^2;
            dist=sqrt(dx+dy);
            [a,b]=min(dist); %index of closest starting node
            count(tic,toc)=count(tic,toc)+1;
            w(tic,toc)=w(tic,toc)+(z(i)-z(b));
            tic=tic+1;
        end
        toc=toc+1; tic=1;
    end
end
%     for j=i+1:length(dat)
%         a(tally,1)=((dat(i,1)-dat(j,1)).^2+(dat(i,2)-dat(j,2)).^2).^.5;
%         a(tally,2)=rad2deg((pi()/2)-atan(((dat(i,2)-dat(j,2))/(dat(i,1)-dat(j,1)))));
%         tally=tally+1;
%     end
% end
snarf=a;

% z=z';
% for i=1:length(nodedat)
%     if i<=length(netdat)
%         a(i,1)=rad2deg((pi()/2)-atan(((nodedat(i,2)-nodedat(netdat(i),2))/(nodedat(i,1)-nodedat(netdat(i),1)))));
%     else
%         a(i,1)=0;
%     end
% end


% a=[];
% tally=1;
% for i=1:length(dat)
%     for j=i+1:length(dat)
%         a(tally,1)=((dat(i,1)-dat(j,1)).^2+(dat(i,2)-dat(j,2)).^2).^.5;
%         a(tally,2)=abs(dat(i,3)-dat(j,3));
%         plot(dat(i,1),dat(i,2),'or',dat(j,1),dat(j,2),'ok'); hold on;
%         a(tally,3)=rad2deg((pi()/2)-atan(((dat(i,2)-dat(j,2))/(dat(i,1)-dat(j,1)))));
%         tally=tally+1;
%     end
% end

        