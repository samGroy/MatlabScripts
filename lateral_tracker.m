function [position,time]=lateral_tracker(basenm,steps)
%lateral_tracker: tracks the lateral migration of a river by finding the
%x position of max(Q) each time step.
position=zeros(steps,1);

filesys='';
filenm= [filesys basenm '.nodes' ];
nfid=fopen(filenm,'r');

filenm=[filesys basenm '.q'];
qfid=fopen(filenm,'r');

for i=1:steps
    tm = fscanf(nfid,'%f',1);
    if i==2, t=tm; end
    tm = fscanf(qfid,'%f',1);
    nn = fscanf(nfid,'%d',1);
    nn = fscanf(qfid,'%d',1);
    n=fscanf(nfid,'%f',[4,nn]);
    q=fscanf(qfid,'%f',[1,nn]);
    
    [val, ind]=max(q);
    position(i)=n(1,ind);
end
time=t*(0:1:steps-1);
fclose(nfid);
fclose(qfid);