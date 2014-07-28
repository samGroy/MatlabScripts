%---convert weak zones produced in the FLAC strain softening model to
%layers in CHILD---%
%--written by Sam Roy and Greg Tucker 6/28/2012--%

%don't run this after you have completed the final CHILD run, bad things happen to the .z file!

%order should be 1) design basic zx plane strength fabric in Surfer, 
%2) use the 'extrudey' script to ready the matrix format (xyzcsirfaceid),
%3) make sure z coordinates are correct and that the model extruded into y,
%4) run CHILD initially for one time step to create surface, 5) run this script, 
%making sure basenm, cname, etc match the name of the flac and CHILD files

%After this script writes the lay file, you can proceed to run CHILD. After
%the run completes, you can view erodibility by running 'creadxyzb' and
%'readlaykbx' in that order, again making sure that basenm matches your
%run name in 'readlaykbx'.

%run empwidth to plot the width calculation per node. Remember to change
%basenm to your respective run name.

clear
basenm='dip50_3';
ts=0;
%read synthesized volume. 
A=importdata('dip50_3.txt',' ');
x=[A(:,1)];
y=[A(:,2)];
z=[A(:,3)];
c=[A(:,4)];
sid=A(:,5);
bottom_of_model = min(z);
top_of_model=max(z);
%this next part is now obsolete if you read in surface data with a surface
%id (sid)
%read synthesized surface, these nodes are included in the volume above.
%B=importdata('testsurf.txt');
%xsurf=[B(:,1)];
%ysurf=[B(:,2)];
%zsurf=[B(:,3)];
%csurf=[B(:,4)];

%this writes the surface coordinates if there is a sid value. 
%A second surface txt file is no longer necessary
for i=1:length(x) 
    
    if sid(i)==1
        xsurf(i,1)=x(i);
        ysurf(i,1)=y(i);
        zsurf(i,1)=z(i);
        csurf(i,1)=c(i);
    end
end

% Read child prelim file. 
cname = 'dip50_3';
xyzb = creadxyzb(cname,1);
% xc = xyzb(:,1);
% yc = xyzb(:,2);
% zc = xyzb(:,3);
% bc = xyzb(:,4);


%now using TriScatteredInterp to assign child node elevation based on "flac" elevation:
F=TriScatteredInterp(xsurf,ysurf,zsurf);
zci = F((xyzb(:,1)),(xyzb(:,2)));
%zci(isnan(zci))=0;
ind=1:length(zci);

%interpolation is required to remove nans. A problem arises for some runs,
%comment out if necessary.
%zci=interp1(ind(~isnan(zci)),zci(~isnan(zci)), ind, 'linear'); 
%ti=0:.1:1;
%[qx,qy]=meshgrid(ti,ti);
%qz=F(qx,qy);
%mesh(qx,qy,qz); hold on;

for i=1:length(xyzb) %clean boundaries to export surface
    if xyzb(i,4)~=0
        zci(i)=0;
    end
end

%comented out for now
% % write zci to .z file
% filenm=[basenm '.z'];
% fid = fopen(filenm,'w');
% fprintf(fid,' %.2f\n',ts);
% fprintf(fid,'%d\n',length(zci));
% for i=1:length(zci)
%     fprintf(fid,'%f\n',zci(i));
% end
% fclose(fid);

for i=1:length(xyzb)%remove boundary nodes!
    if xyzb(i,4)==0
        zcii(i)=zci(i);
        xc(i) = xyzb(i,1);
        yc(i) = xyzb(i,2);
        zc(i) = xyzb(i,3);
        bc(i) = xyzb(i,4);
    end
end

% for i=1:length(xyzb) %remove boundary nodes!
%     if xyzb(i,4)==0
%         xc(i) = xyzb(i,1);
%         yc(i) = xyzb(i,2);
%         zc(i) = xyzb(i,3);
%         bc(i) = xyzb(i,4);
%     end
% end

% plot flac and child nodes, plot every 10th child node if rendering is slow
plot3(x,y,z,'+');
hold on
plot3(xc,yc,zcii,'.r');
hold on



% DEFINE AND INITIALIZE VARIABLES NEEDED FOR LAYERS
dz=200; %vertical step taken to find layer transitions
num_child_nodes = length(xc);
max_num_layers = 3;   % maximmum number of layers at a node
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
    for z_current = zcii(i):-dz:bottom_of_model
        zcounter=zcounter+1;
        zfromsurf=z_current-zcii(i);
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
            baselay(nlayers(1,i)-1,i)=zcii(i)+zcdiff(n-1,num_child_nodes+1); %find depth just above the transition, then subtract it from the respective child elevation (zci) to get the base-of-layer elevation
            layer_cohesion(nlayers(1,i)-1,i)=zcmatrix(n,i); %assign cohesion of layer
        end
    end
    layer_thickness(1,i)=zcii(i)-baselay(1,i); %assigns top layer thickness
    layer_thickness(2:end,i)=baselay((2:end)-1,i)-baselay(2:end,i); %assigns thicknesses to all layers
end
nlayers=nlayers-1; %return to correct sum of layers

%convert cohesion to erodibility
alpha = 100;
epsilon = 1;
kb = alpha ./ ( layer_cohesion + epsilon );

% % make baselay visible as points plot every 10 points if rendering is slow 
% for i=1:length(xc)
%     for j=1:nlayers(i)
%         plot3(xc(i),yc(i),baselay(j,i),'.g');
%         hold on;
%     end
% end
% hold off;

%attempt to interpolate surfaces from baselay, faster rendering compared to
%point data
giggity=baselay';
coh=nan(length(giggity),max_num_layers);
xx=xc'; yy=yc';
snarf=zeros(length(giggity),1);
laysurf=nan(length(giggity),max_num_layers);
for i=1:length(giggity)
    snarf(i)=sum(~isnan(giggity(i,:)));
    if snarf(i)==1
        laysurf(i,max_num_layers)=giggity(i,1);
        coh(i,3)=layer_cohesion(1,i);
    end
    if snarf(i)==2
        laysurf(i,2)=giggity(i,1);
        laysurf(i,3)=giggity(i,2);
        coh(i,2)=layer_cohesion(1,i);
        coh(i,3)=layer_cohesion(2,i);
    end
    if snarf(i)==3
        laysurf(i,1)=giggity(i,1);
        laysurf(i,2)=giggity(i,2);
        laysurf(i,3)=giggity(i,3);
        coh(i,1)=layer_cohesion(1,i);
        coh(i,2)=layer_cohesion(2,i);
        coh(i,3)=layer_cohesion(3,i);
    end
end

tri=delaunay(xx,yy);
trisurf(tri,xx,yy,laysurf(:,1),coh(:,1));
hold on;
trisurf(tri,xx,yy,laysurf(:,2),coh(:,2));
hold on;
trisurf(tri,xx,yy,laysurf(:,3),coh(:,3));
hold off;



%--now, make the lay file--%
% Create layer file and open it for writing (first make sure it doesn't
% already exist)


filenm= [ basenm '.lay' ];
lfid=fopen(filenm,'r');
% if lfid>0
%     fprintf('Layer file "%s" already exists.\n',filenm);
%     error('Try a different name.\n');
% end
lfid = fopen(filenm,'w');
if lfid<1
    error('Unable to create layer file.\n');
end

fprintf('now writing lay file.\n');

ALLUVIUM_FLAG=1;

% Write header information
fprintf(lfid,' %.2f\n',ts);
fprintf(lfid,'%d\n',length(bc));

% Loop over nodes to write information
for j=1:length(xc)
    
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',nlayers(j));

    % For each layer at this node, write information
    for i=1:nlayers(j)
        
        % Basic properties: creation time, recent activity time, exposure
        % time, thickness, erodibility, and regolith/bedrock flag
        fprintf(lfid,'%.2f %.2f %.2f\n',0,0,0);
        fprintf(lfid,'%.2f %f %.0f\n',layer_thickness(i,j),kb(i,j),0);
        fprintf(lfid,'%.2f\n',layer_thickness(i,j));
        
%         % Grain size information
%         if numg==1
%             fprintf(lfid,'%.2f\n',laydat(j,i,1));   % this line gives 'old' format: if one size, just write thickness
%         elseif numg>1
%             fprintf(lfid,'%.2f ',laydat(j,i,1)-sum(laydat(j,i,5:(3+numg))));  % this line gives 'old' format: if more than one, write thickness of size 1 as total thick minus sum of sizes 2-N
%             for k=2:numg-1
%                 fprintf(lfid,'%f ',laydat(j,i,4+(k-1)) );
%             end
%             fprintf(lfid,'%f\n',laydat(j,i,4+(numg-1)) );
%         end
        
    end
    
end
fprintf('done.\n');