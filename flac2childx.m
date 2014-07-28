%flac2child.m: reads FLAC output and converts to:
% 1 - CHILD uplift rate map
% 2 - file with erodibility info

% GT Sep '11

% Read flac file
fname = 'voxler_zone_88x.txt';
fid = fopen(fname,'r');
nflacnodes = fscanf(fid,'%d',1);
fprintf('There are %d nodes in flac run\n',nflacnodes);
flacdata = fscanf(fid,'%f',[6,nflacnodes]);
fclose(fid);
id = flacdata(1,:);
xf = flacdata(2,:);
yf = flacdata(3,:);
zf = flacdata(4,:);
uf = flacdata(5,:);
cf = flacdata(6,:);

% Read child prelim file
cname = 'ChildRun/testf2c1_setup';
xyzb = creadxyzb(cname,1);
xc = xyzb(:,1);
yc = xyzb(:,2);
zc = xyzb(:,3);
bc = xyzb(:,4);

% Generate a set of tiny random offsets for the flac nodes. This is needed
% to interpolate from flac to child -- if the grid is perfectly regular
% (rectilinear), the interpolation algorithm chokes.
tiny_offset_x = 0.001*rand(size(xf));
tiny_offset_y = 0.001*rand(size(yf));

% interpolate uplift rates from flac to child, and write to file
uc=griddata(xf+tiny_offset_x,yf+tiny_offset_y,uf,xc,yc);
uc( find( isnan( uc ) ) ) = 0.0;
fid = fopen('ChildRun/flac_uplift_map001','w');
for i=1:length(uc)
    fprintf(fid,'%f\n',uc(i));
end
fclose(fid);

% interpolate cohesion and turn it into erodibility
cc=griddata(xf+tiny_offset_x,yf+tiny_offset_y,cf,xc,yc,'nearest');
alpha = 100000;   % Converts c in Pa to Kb
epsilon = 1;   % Minimum c in Pa (prevents division by zero)
default_cohesion = 5e7;
default_kb = alpha ./ default_cohesion;
cc( find( isnan( cc ) ) ) = default_cohesion;
cerody = alpha ./ ( cc + epsilon );

% write erodibility to file
fid = fopen('ChildRun/flac_erodibility_map001','w');
for i=1:length(cerody)
    fprintf(fid,'%.7f\n',cerody(i));
end
fclose(fid);





