%read synthesized volume. This synthesized volume has 15 nodes.
A=importdata('test.txt',' ');
x=[A(:,1)];
y=[A(:,2)];
z=[A(:,3)];
c=[A(:,4)];
bottom_of_model = min(z);

%read synthesized surface, these nodes are included in the volume above.
%There ar 5 for each vector.
B=importdata('testsurf.txt');
xsurf=[B(:,1)];
ysurf=[B(:,2)];
zsurf=[B(:,3)];
csurf=[B(:,4)];

%parameters used for the interpolation below. I think the step size (nl)
%should equal the node spacing for the CHILD grid?
% ==> I DON'T THINK YOU NEED xl OR yl
cstrong=max(c);
%xl=.1;
%yl=.1;
dz=.1;

%---I'm forgoing this step for now, just using a mesh created by
%TriScatteredInterp---
% Read child prelim file. A problem occurs: if I run this script twice, bc
% is taken out of the matrix for some unknown reason (says index is out of bounds).
% ==> FOR NOW, WE'LL JUST MAKE UP A SET OF CHILD POINTS
xc=[0.2 0.51 0.8];
yc=[0.31 0.49 0.76];
zc=[0 .01 .02];
bc=[1 0 1];
num_child_nodes = length(xc);
%cname = 'test01';
%xyzb = creadxyzb(cname,1);
%xc = xyzb(:,1);
%yc = xyzb(:,2);
%zc = xyzb(:,3);
%bc = xyzb(:,4);

% Generate a set of tiny random offsets for the flac nodes. This is needed
% to interpolate from flac to child -- if the grid is perfectly regular
% (rectilinear), the interpolation algorithm chokes.
% ==> WE MIGHT NOT NEED THIS ANYMORE IF TRISCATTEREDINTERP HANDLES
% REGULARLY SPACED DATA
tiny_offset_x = 0.001*rand(size(x));
tiny_offset_y = 0.001*rand(size(y));

%---I'm forgoing this step for now---
% interpolate z from flac to child, and write to .txt file
%zinterp=griddata(x+tiny_offset_x,y+tiny_offset_y,z,xc,yc);
%zinterp( find( isnan( zinterp ) ) ) = 0.0;
%fid = fopen('test01z.txt','w');
%for i=1:length(zinterp)
    %fprintf(fid,'%f\n',zinterp(i));
%end
%fclose(fid);

%now using TriScatteredInterp:
F=TriScatteredInterp(xsurf,ysurf,zsurf);
zci = F(xc,yc);
%ti=0:.1:1;
%[qx,qy]=meshgrid(ti,ti);
%qz=F(qx,qy);
%mesh(qx,qy,qz); hold on;
plot3(x,y,z,'+');
hold on
plot3(xc,yc,zci,'.')
hold off

%find euclidean distance from probe distance xl*n to synthesized node,
%determine the closest flac node and use that node's cohesion value (not there yet).

% DEFINE AND INITIALIZE VARIABLES NEEDED FOR LAYERS
max_num_layers = 3;   % maximmum number of layers at a node
nlayers = ones(1,num_child_nodes);  % number of layers at each child node
layer_thickness = nan(max_num_layers,num_child_nodes);  % thickness of each layer
layer_cohesion = nan(max_num_layers,num_child_nodes);  % cohesion of each layer
layer_cohesion(1,:)=cstrong; %I put this in to fill the first row, but it assumes that no weak zones intersect the surface

% FOR EACH CHILD NODE (x,y) LOCATION
for i=1:num_child_nodes
   
    % FOR EACH DEPTH INCREMENT
    for z_current = zci:-dz:bottom_of_model
        
        % FIND CLOSEST FLAC NODE TO THE CURRENT LOCATION
        distsq=((xc(i)-x).^2+(yc(i)-y).^2+(z_current-z).^2);
        [mindistsq,closest_flac_node_id] = min(distsq);
        fprintf('The closest flac node to %f,%f,%f is node %d, whose coords are %f,%f,%f\n',...
            xc(i),yc(i),z_current,closest_flac_node_id,x(closest_flac_node_id),...
            y(closest_flac_node_id),z(closest_flac_node_id) );
        
        % layer cohesion assignment algorithm
        if c(closest_flac_node_id)~=cstrong
            %we need to base this on a transition, otherwise, every
            %z_current probe will create a new layer, not what we want...
            nlayers(:,i)=nlayers(:,i)+1;
            layer_cohesion(nlayers(:,i),i)=c(closest_flac_node_id);
            baselay(nlayers(:,i)-1)=z_current;
        end
    end
    
end

% for zi=max(z):zl:min(z)
%     for xi=min(x):xl:max(x)
%         for yi=min(y):yl:max(y)
%             distsq=((xi-x).^2+(yi-y).^2+(zi-z).^2);
%             [mindist,loc]=min(distsq);
%         end
%     end
% end