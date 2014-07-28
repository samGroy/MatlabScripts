function [t,v] = cvolplotuniform( fname, n,tm, total_area )
%CVOLPLOT: Plots volume above base level versus time from a Child run. Run must
%          have used OPTTSOUTPUT option.
%    Usage: [t,v] = cvolplot( fname, n )
%      fname = file name (including path)
%      n = number of lines in <fname>.storm and <fname>.vols files
%    Returns:
%      t = time vector (x-axis)
%      v = volume vector (y-axis)
%   GT, May 2002
storm=1;
volumefile = [fname '.vols'];
vfid = fopen( volumefile, 'r' );
if vfid <= 0, error(['Unable to open ' volumefile]);end
stormfile = [fname '.storm'];
sfid = fopen( stormfile, 'r' );
if sfid <= 0, storm=0; end
ufile = [fname '.up'];
ufid = fopen( ufile, 'r' );
if ufid <= 0, error(['Unable to open ' ufile]);end

% If user gives zero or negative for n, figure it out from file
if nargin<2 || n<=0
   mycomputer = computer;
   if strncmp(mycomputer,'PCWIN',5)
       error('Cannot automatically determine file length on a windows machine.');
   end
   [s,w]=unix(['wc -l ' fname '.vols']);
   n = str2num( w(1:9) );
end

% If user didn't specify area, read it from file
varea = cread([fname '.varea'],1);  % assumes area has stayed constant!
total_area = sum(varea);



v = fscanf( vfid, '%f', [1,n] );
fclose( vfid );
if storm>0
    storm=(fscanf(sfid, '%f', [3,n]))';
    t=storm(:,1)+storm(:,3);
    t=cumsum(t);
else
    iteration=tm/n;
    olditeration=iteration;
    t=nan(n,1);
    for i=1:n;
        iteration=olditeration*i;
        t(i)=iteration;
    end
end
tm=fscanf(ufid, '%f',1);
nn=fscanf(ufid, '%f',1);
junk=fscanf(ufid,'%f',nn+2);
u=fscanf(ufid, '%f',1);

% Sometimes the two files won't have equal length, because of output
% buffering or copying of files during an active run. In this case, set the
% length to the shorter of the two (if v is shorter, this is already taken
% care of, so we need only handle the case in which t is shorter than v).
if length(t)<n
    n=length(t);
    v=v(1:n);
end
vt=((t*u+100)*total_area)-(v');
% vt=(1000*total_area)-(v');
% Divide volume by total area to get mean altitude and plot
%  melev=v./total_area;
plot( t, v )
grid on

