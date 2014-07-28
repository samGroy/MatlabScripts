function m=tortuositymap(basenm,ts)
%basenm='fractureset05';
filesys=[''];
%ts=5;
filenm= [filesys basenm '.net' ];
netfid=fopen(filenm,'r');
if netfid==0, error(['Unable to open ' filenm]);end
filenm= [filesys basenm '.nodes' ];
nodfid=fopen(filenm,'r');
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
end

nodedat=nodedat';
netdat=(netdat+1)';
trig=nan(intnodes,4);
tort=zeros(allnodes,3);

for ii=1:intnodes
   if nodedat(ii,1)<nodedat(netdat(ii),1) && nodedat(ii,2)<nodedat(netdat(ii),2) % quadrant I
       trig(ii,2)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))));
   end
   if nodedat(ii,1)>nodedat(netdat(ii),1) && nodedat(ii,2)<nodedat(netdat(ii),2) % quadrant II
       trig(ii,2)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))))+pi/2;
   end
   if nodedat(ii,1)>nodedat(netdat(ii),1) && nodedat(ii,2)>nodedat(netdat(ii),2) % quadrant III
       trig(ii,2)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))))+2*pi/2;
   end
   if nodedat(ii,1)<nodedat(netdat(ii),1) && nodedat(ii,2)>nodedat(netdat(ii),2) % quadrant IV
       trig(ii,2)=abs(atan((nodedat(netdat(ii),2)-nodedat(ii,2))/(nodedat(netdat(ii),1)-nodedat(ii,1))))+3*pi/2;
   end
   trig(ii,1)=((nodedat(netdat(ii),2)-nodedat(ii,2)).^2+(nodedat(netdat(ii),1)-nodedat(ii,1)).^2).^.5;
   trig(ii,3)=netdat(ii);
   if trig(ii,3)>intnodes
       trig(ii,4)=trig(ii,3);
   else
   trig(ii,4)=netdat(netdat(ii));
   end
end

for iii=1:intnodes
    if netdat(iii)>intnodes
        tort(iii,1)=tort(iii,1);
    else
    tort(iii,1)=trig(iii,2)-trig(netdat(iii),2);
    end
    tort(iii,2)=((nodedat(trig(iii,4),2)-nodedat(iii,2)).^2+(nodedat(trig(iii,4),1)-nodedat(iii,1)).^2).^.5;
    tort(iii,3)=abs(tort(iii,1)/tort(iii,2));
end
m=tort(:,3);