function a=waveorient1(basenm,ts)
filesys=[''];
%ts=5;
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.nodes' ];
nodfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
filenm= [filesys basenm '.tri' ];
tfid=fopen(filenm,'r');

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
  tm=fscanf(zfid,'%f',1);
  allnodes= fscanf(zfid,'%d',1);
  z=fscanf(zfid,'%f',[1,allnodes]);
  nt = fscanf(tfid,'%d',1);
  t=fscanf(tfid,'%f',[9,nt]); 

end

tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
nodedat=nodedat';
netdat=(netdat+1)';
z=z';
for i=1:length(nodedat)
    if i<=length(netdat)
        a(i,1)=rad2deg((pi()/2)-atan(((nodedat(i,2)-nodedat(netdat(i),2))/(nodedat(i,1)-nodedat(netdat(i),1)))));
    else
        a(i,1)=NaN;
    end
end
v

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

        