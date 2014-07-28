function [x,y,s,b,z]=creadxyzb2(basenm,ts)
% CREADXYZ: Reads node x,y,z coordinates and boundary flag ("b")
%           for one timestep of a CHILD run
%           and returns them in a Nx4 matrix, where N=number of nodes.
%
%  Usage: m = creadxyzb( basenm, ts )
%
%  Parameters:
%    filenm -- name of edge file
%    ts -- time slice # to fetch data for
%
filesys=[''];
filenm= [filesys basenm '.nodes' ];
nfid=fopen(filenm,'r');
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
for i=1:ts
  tm = fscanf(nfid,'%f',1);
  nn = fscanf(nfid,'%d',1);
  n=fscanf(nfid,'%f',[4,nn]);
  tm = fscanf(zfid,'%f',1);
  nn = fscanf(zfid,'%d',1); 
  z=fscanf(zfid,'%f',[1,nn]);
end
x = rot90(n(1,:),3);
y = rot90(n(2,:),3); 
s = rot90(n(3,:),3); 
b = rot90(n(4,:),3);
z = rot90(z,3); 
fclose(nfid);
fclose(zfid);


