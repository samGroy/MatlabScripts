function [nodedat,h]=tortuositymap2(basenm,ts, nwavelength)

filesys=[''];
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.nodes' ];
nodfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
if zfid<=0, error('Unable to open elevation file'),end
filenm= [filesys basenm '.tri' ];
tfid=fopen(filenm,'r');
if tfid<=0, error('Unable to open triangle file'),end
%nd(:,3)=nd(:,3)+5;
fprintf('CPLOTAREA: Reading data ...\n');

for i=1:ts
  tm=fscanf(netfid,'%f',1);
  fprintf('Time slice %d (T=%f)\n',i,tm);
  tm=fscanf(nodfid,'%f',1);
  tm = fscanf(tfid,'%f',1);
  intnodes= fscanf(netfid,'%d',1);
  nt = fscanf(tfid,'%d',1);
  t=fscanf(tfid,'%f',[9,nt]); 
  allnodes= fscanf(nodfid,'%d',1); 
  netdat=fscanf(netfid,'%d',[1,intnodes]);
  nodedat=fscanf(nodfid,'%f',[4,allnodes]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,allnodes]); 
end

nodedat=nodedat';
nodedat(:,5)=0;
netdat=(netdat+1)';
steparc=nan(allnodes,3);
arclength=nan(nwavelength,3);
tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
z=z';

for i=1:length(netdat);
    a=1; b=i;
    while a<=nwavelength
        if nodedat((netdat(b)),4)==2, break, end
        dx=nodedat(b,1)-nodedat(netdat(b),1); dx=dx.^2;
        dy=nodedat(b,2)-nodedat(netdat(b),2); dy=dy.^2;
        arclength(a,1)= nodedat(b,1);
        arclength(a,2)= nodedat(b,2);
        arclength(a,3)=sqrt(dx+dy);
        if a==round(nwavelength/2), bb=b; end
        if a==nwavelength, break, end
        a=a+1; b=netdat(b);
    end
    steparc(i,1)=arclength(round(nwavelength/2),1);
    steparc(i,2)=arclength(round(nwavelength/2),2);
    dx=arclength(1,1)-arclength(nwavelength,1); dx=dx.^2;
    dy=arclength(1,2)-arclength(nwavelength,2); dy=dy.^2;
    nodedat(bb,5)=1-(sqrt(dx+dy)/sum(arclength(:,3)));
end

trisurf(tri,nodedat(:,1),nodedat(:,2),z,nodedat(:,5));
colormap jet;
shading flat;
caxis([0 1]);
colorbar;
axis([min(nodedat(:,1)) max(nodedat(:,1)) min(nodedat(:,2)) max(nodedat(:,2)) min(z) max(z)+1000]);