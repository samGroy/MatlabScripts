%---convert weak zones produced in the FLAC strain softening model to
%layers in CHILD---%
%--written by Sam Roy and Greg Tucker 6/28/2012--%

%read synthesized volume. This synthesized volume has 15 nodes.
A=importdata('test.txt',' ');
x=[A(:,1)];
y=[A(:,2)];
z=[A(:,3)];
c=[A(:,4)];
sid=A(:,5);
bottom_of_model = min(z);
top_of_model=max(z);
%read synthesized surface, these nodes are included in the volume above.
%There ar 5 for each vector.
%B=importdata('testsurf.txt');
%xsurf=[B(:,1)];
%ysurf=[B(:,2)];
%zsurf=[B(:,3)];
%csurf=[B(:,4)];
for i=1:length(x) %this writes the surface coordinates if there is a sid value. A second surface txt file is no longer necessary
    if sid(i)==1
        xsurf(i,1)=x(i);
        ysurf(i,1)=y(i);
        zsurf(i,1)=z(i);
        csurf(i,1)=c(i);
    end
end

% Read child prelim file. 
cname = 'test01';
xyzb = creadxyzb(cname,1);
xc = xyzb(:,1);
yc = xyzb(:,2);
zc = xyzb(:,3);
bc = xyzb(:,4);
%xc(find(xc<0))=0;
%yc(find(yc<0))=0;
num_child_nodes = length(xc);

%now using TriScatteredInterp to assign child node elevation based on "flac" elevation:
F=TriScatteredInterp(xsurf,ysurf,zsurf);
zci = F(xc,yc);
%zci(isnan(zci))=0;
ind=1:length(zci);
zci=interp1(ind(~isnan(zci)),zci(~isnan(zci)), ind, 'linear'); %interpolation is required to remove nans
%ti=0:.1:1;
%[qx,qy]=meshgrid(ti,ti);
%qz=F(qx,qy);
%mesh(qx,qy,qz); hold on;
plot3(x,y,z,'+');
hold on
plot3(xc,yc,zci,'.')
hold off

%---I'm forgoing this step for now---
% write zci to .z file
%fid = fopen('test01.z','w');
%for i=1:length(zinterp)
    %fprintf(fid,'%f\n',zinterp(i));
%end
%fclose(fid);

% DEFINE AND INITIALIZE VARIABLES NEEDED FOR LAYERS
dz=.1; %vertical step taken to find layer transitions
max_num_layers = 4;   % maximmum number of layers at a node
nlayers = ones(1,num_child_nodes);  % number of layers at each child node
layer_thickness = nan(max_num_layers,num_child_nodes);  % thickness of each layer
layer_cohesion = nan(max_num_layers,num_child_nodes);  % cohesion of each layer
%layer_cohesion(1,:)=cstrong; %I put this in to fill the first row, but it assumes that no weak zones intersect the surface
baselay=nan(max_num_layers,num_child_nodes); %this keeps track of the depth to each layer
zcmatrix=zeros(top_of_model-bottom_of_model/dz+1,num_child_nodes+1); %use this to keep track of cohesion values with respect to depth-stepping
zcounter=0; %this is used to assign the row in which cohesion and depth will be recorded.
zcdiff=zeros((top_of_model-bottom_of_model/dz)+1,num_child_nodes+1);
% FOR EACH CHILD NODE (x,y) LOCATION
for i=1:num_child_nodes
   zcounter=0;
   
    % setting up the matrix of c and z data used to find layer transitions
    for z_current = zci(i):-dz:bottom_of_model
        zcounter=zcounter+1;
        zfromsurf=z_current-zci(i);
        zcmatrix(zcounter,num_child_nodes+1)=zfromsurf;
        ztester(zcounter,i)=z_current;
        % FIND CLOSEST FLAC NODE TO THE CURRENT LOCATION
        distsq=((xc(i)-x).^2+(yc(i)-y).^2+(z_current-z).^2);
        [mindistsq,closest_flac_node_id] = min(distsq);
        fprintf('The closest flac node to %f,%f,%f is node %d, whose coords are %f,%f,%f\n',...
            xc(i),yc(i),z_current,closest_flac_node_id,x(closest_flac_node_id),...
            y(closest_flac_node_id),z(closest_flac_node_id) );
        zcmatrix(zcounter,i)=c(closest_flac_node_id); %this assigns cohesion per depth for node i

    end

end
zcmatrix(((top_of_model-bottom_of_model)/dz)+2,1:num_child_nodes)=0; %terminate the base of the volume with an exorbitantly large cohesion difference
zcmatrix(((top_of_model-bottom_of_model)/dz)+2,num_child_nodes+1)=-top_of_model+bottom_of_model-dz; %z depth just below the bottom of the model to allow layer termination

%use changes in cohesion from zcmatrix to determine layer changes
for i=1:num_child_nodes
    for n=1:((top_of_model-bottom_of_model)/dz)+1 %can maybe replace this outer loop with diff function
        zcdiff(n,i)=abs(zcmatrix(n,i)-zcmatrix(n+1,i)); %populate zcdiff with cohesion differences
        zcdiff(n,num_child_nodes+1)=zcmatrix(n+1,num_child_nodes+1); %add depth on furthest column
        
        if zcdiff(n,i)>=1 %where a cohesion transition greater than 1 occurs
            nlayers(1,i)=nlayers(1,i)+1; %add a layer
            baselay(nlayers(1,i)-1,i)=zci(i)+zcdiff(n-1,num_child_nodes+1); %find depth just above the transition, then subtract it from the respective child elevation (zci) to get the base-of-layer elevation
            layer_cohesion(nlayers(1,i)-1,i)=zcmatrix(n,i); %assign cohesion of layer
        end
    end
    layer_thickness(1,i)=zci(i)-baselay(1,i); %assigns top layer thickness
    layer_thickness(2:end,i)=baselay((2:end)-1,i)-baselay(2:end,i); %assigns thicknesses to all layers
end
nlayers=nlayers-1;
