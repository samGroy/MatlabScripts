function [id,x,y,z,v,c,er,ps] = readflac2( filenm )
%READFLAC2: Reads in data from a FLAC run.
%  Usage: [id,x,y,z,v,c] = readflac2( filenm )
%  Inputs:  filenm - name of FLAC file
%  Returns: 6 x N matrix containing ID number, x, y, z position,
%  vertical velocity, and cohesion.
%
%  The assumed file format has the number of nodes on the first line, and
%  then on each subsequent line, separated by spaces, the ID, x, y, z, 
%  v, and c values for each node.
%
% Version 1
% GT, for Taiwan project, Oct 2007; Nov 2010
% edited by SGR for NZ proposal 12/19/12
% Open the file
fid = fopen( filenm, 'r' );
if fid<=0
    error(['Unable to open FLAC file "' filenm '"']);
end

% Read the number of nodes
NN = fscanf(fid,'%d\n',[1]);

% Read the rest of the data
m = fscanf(fid,'%f',[8,NN]);
id=m(1,:);
x=m(2,:);
y=m(3,:);
z=m(4,:);
v=m(5,:);
c=m(6,:);
er=m(7,:);
ps=m(8,:);

% Close the file
fclose(fid);

