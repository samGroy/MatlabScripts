function product = synthlay_grad5( basenm, matrix, ts, step, numg, gfrac1, gfrac2, gfrac3, gfrac4)
%SGR 7/2012
%used to create MC fracture sets, intersecting faults produced by strain
%that can be eroded preferentially by surface processes.
%faults can cross along strike, need to specify in file, origin
%coordinates, dip, and strike for two planes.

% SGR 9/2012 UPDATE
% Layer thickness algorithm made more efficient and capable of handling
% many more faults, front end allows for 6 faults now, but more could be
% added if desired. All faults are uniform single cohesion.

% SGR 10/2012 UPDATE
% Added input functions for grain size info.

% SGR 5/2013 THOUGHT 
% why not indicate number of faults, then run a for loop instead of having
% to define 11 faults each time? Contain fault data in a single matrix. 
% This would shorten the code length and add flexibility to the number of 
% faults that can be generated. The fault matrix needs x0, y0, z0, strike,
% dip, thickness. Cohesion values can't vary between faults.

% SGR 5/2013 ANTECEDENCE alternative
% See Synthlay_grad4 for antecedence algo.

% SGR 7/2013 Pre-UPDATE
% Plan to carry all fault data in 3D matrix. Column: variables, row:
% faults, 3rd D: number of runs. Did this to make code more compact, adding
% loops and reducing the massive number of variables. Hopefully easier to
% read, may be faster? Also allows for however many faults you want.
% Convert degree input to rad value
% basenm='asymflt01_1'; step=1;
% ts=0;
% step=1;
% numg=1;

geometry=matrix(2:length(matrix((matrix(:,1)~=0),1)),:);
cohesion=matrix(1,1:4);

maxfaults=length(geometry(:,1));
maxplanes=maxfaults*6+1;

dipangle=degtorad(geometry(:,1));
strikeangle=degtorad(geometry(:,2));

% Calculate vertical thickness of fault plane
w_mult=5; %multiplier determines total length of fault core and damage zone, input 'fault_thickness*' is fault core width
fthick=abs((geometry(:,3)*w_mult)/cos(dipangle)); %vertical thickness of fault plane


% p0 is the origin point defined from the geometry matrix, xyz data cols 5-7
% p1 is the point taken down strike from the origin
% p2 is the point taken down dip from the origin
% p0=geometry(:,5:7);
% p1=[geometry(:,5)+cos(strikeangle),geometry(:,6)+sin(strikeangle),geometry(:,7)]; 
% p2=[geometry(:,5)+cos(strikeangle+(pi/2)),geometry(:,6)+sin(strikeangle+(pi/2)),slope+geometry(:,7)]; 

% the vector calcs are just p1-p0 and p2-p0, the ortho calc is the cross
% product between vector 1 and 2
vector1=[geometry(:,5)+cos(strikeangle),geometry(:,6)+sin(strikeangle),geometry(:,7)]-geometry(:,5:7);
vector2=[geometry(:,5)+cos(strikeangle+(pi/2)),geometry(:,6)+sin(strikeangle+(pi/2)),tan(dipangle)+geometry(:,7)]-geometry(:,5:7);
ortho=cross(vector1,vector2);

% Read child prelim file, omit boundary nodes
cname = basenm;
xyzb = creadxyzb(cname,step);
xc = rot90(xyzb((xyzb(:,4)==0),1));
yc = rot90(xyzb((xyzb(:,4)==0),2));
zc = rot90(xyzb((xyzb(:,4)==0),3));

% For loop to define the elevation of all points defining all planes. There
% are 6 planes per fault currently. The number of nodes equals the number of interior nodes arranged by CHILD in the above step. There is an option to make fault zones
% symmetric or asymmetric, currently asymmetric faults weaken to the
% hanging wall (the bottom plane).
depth_to_planes=zeros(maxplanes+1,length(xc)); % keeps track of depth to every plane at any point
plane_id=depth_to_planes*nan; 
% plane_id allows for fault cores to supercede other layers with greater cohesion,
% therefore at any given point in 3D, if weak and strong layers overprint,
% the weak layer prevails. plane_id also keeps track of each sublayer in each
% fault, so as not to lose track of the sublayer order, especially when
% several faults overlap at a point.
depth_to_planes(1,:)=zc; 
depth_to_planes(maxplanes+1,:)=-1000000; 
plane_id(1,:)=5;
plane_id(maxplanes+1,:)=0;
jold=0;
for i=1:maxfaults %bdfeca,024531, add 1, 4/5, 3/5, 2/5, 1/5, 0 to the original plane
    plane_calc=((ortho(i,1)*(xc-geometry(i,5))+ortho(i,2)*(yc-geometry(i,6)))/ortho(i,3))+geometry(i,7); %the plane calc is the equation for a 3D plane.
    if geometry(i,4)==0 % symmetric fault zones
        for j=1:6
            depth_to_planes(j+jold+1,:)=plane_calc+fthick(i,1)*(1-((j-1)/5)); % vertically space the 6 planes defining the fault zone.
            depth_to_planes(j+jold+1,(depth_to_planes(j+jold+1,:)>zc))=zc((depth_to_planes(j+jold+1,:)>zc)); %omit supersurface plane nodes.
            plane_id(j+jold+1,:)= j-1;
        end
    elseif geometry(i,4)==1 % asymmetric fault zones
        for j=1:6
            if j>3
                depth_to_planes(j+jold+1,:)=plane_calc;
            else
                depth_to_planes(j+jold+1,:)=plane_calc+fthick(i,1)*(1-((j-1)/3));
            end
            depth_to_planes(j+jold+1,(depth_to_planes(j+jold+1,:)>zc))=zc((depth_to_planes(j+jold+1,:)>zc)); %omit supersurface plane nodes.
            plane_id(j+jold+1,:)= j-1;
        end
    end
    jold=j+jold;
end

% Introduce empty matrices that will hold data for sublayer thickness,
% cohesion, number of layers at each node. 
thickness=zeros(maxplanes, length(xc));
C=nan(maxplanes,length(xc));
nlayers=zeros(1,length(xc));

% Sort depth_to_planes from highest to lowest, and use the index to sort
% plane_id so we know which plane exists where, which overlap, etc.
[Y,I]=sort(depth_to_planes,1,'descend');
% col=1:1:length(xc);
depth_to_planes=Y;
% plane_id=plane_id(I,col); crashes
for i=1:length(xc)
    plane_id(:,i)=plane_id(I(:,i),i);
end
% This algorithm keeps track of all of the sublayers for each fault and
% host rock. Thicknesses of each sublayer are measured and put in the
% thickness matrix. In the case of overlapping sublayers, the lowest
% cohesion sublayer overprints all others, and those with equal cohesion are
% consolidated as a single sublayer. Fault cores will therefore never
% be overprinted, and host rock will always be overprinted in the
% presence of faults.
% n=2; tic=1; tally=0; sally=tally; bally=sally; % tally keeps track of the outer layer of the fault between planes a and b, sally keeps track of the internal layers between c and d, bally keeps track of the fault core between e and f
for i=1:length(xc)
    n=2; tic=1; tally=0; sally=tally; bally=sally; % tally keeps track of the outer layer of the fault between planes a and b, sally keeps track of the internal layers between c and d, bally keeps track of the fault core between e and f
    while n <= maxplanes+1 && sum(thickness(:,i))<-depth_to_planes(maxplanes+1,i)+100 % for a depth between the top and bottom of the model domain.
        if plane_id(n,i)==0 %top of outer fault layer, bottom of host layer
            if tally==0 && sally==0 && bally==0 % this signifies the end of a host rock layer.
                thickness(tic,i)=depth_to_planes(n-1,i)-depth_to_planes(n,i); % thickness is depth difference between planes
                C(tic,i)=cohesion(4); % assign host rock cohesion value to sublayer
                tic=tic+1; % move on to next layer, transition to a new layer with a new cohesion value and thickness.
            end
            if tally>0 && sally==0 && bally==0 % Indicates the outermost sublayer of a fault. Accumulate the outer sublayer thicknesses.
                thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
                C(tic,i)=cohesion(3); % outermost fault sublayer cohesion value
            end
            tally=tally+1; % when tally > 0, the sublayer is at least the outermost layer of a fault.
            n=n+1; continue % move on to the next index in order.
        end
        
        if plane_id(n,i)==1 %top of middle fault layer, bottom of outer fault layer
            if sally==0 %terminate the outermost fault sublayer thickness, the lower cohesion middle fault layer supercedes.
                thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
                C(tic,i)=cohesion(3); % assign the outermost fault sublayer cohesion value to the terminated layer
                tic=tic+1; % move on to next sublayer
            end
            if sally>0 && bally==0 % accumulate the middle layer thicknesses
                thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
                C(tic,i)=cohesion(2); % assign middle fault sublayer cohesion value
            end
            sally=sally+1; % When sally > 0, the sublayer is at least the middle layer of the fault
            n=n+1; continue
        end
        
        if plane_id(n,i)==2 %top of core fault layer, bottom of middle fault layer, indicates that fault core is met at this depth, all other layers are superceeded.
            if bally==0 %terminate the middle sublayer fault thickness
                thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
                C(tic,i)=cohesion(2);
                tic=tic+1;
            end
            if bally>0 % accumulate the core layer thicknesses
                thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
                C(tic,i)=cohesion(1); % assign fault core cohesion value
            end
            bally=bally+1; % when bally > 0, the sublyer is at least the fault core.
            n=n+1; continue
        end
        %         layer bottoms
        if plane_id(n,i)==3 %bottom of core layer, top of middle layer
            bally=bally-1; % the bottom side of the fault core
            thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
            C(tic,i)=cohesion(1);
            n=n+1;
            if bally==0 % if all overlapping fault cores are terminated at this depth,
                tic=tic+1; %  move to next sublayer thickness
            end
        end
        if plane_id(n,i)==4 %bottom of middle layer, top of outer layer
            sally=sally-1; % the bottom side of the middle fault sublayer
            thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
            C(tic,i)=cohesion(2);
            n=n+1;
            if sally==0
                tic=tic+1;
            end
        end
        
        if plane_id(n,i)==5 %bottom of outer layer, top of host
            tally=tally-1;
            thickness(tic,i)=thickness(tic,i)+depth_to_planes(n-1,i)-depth_to_planes(n,i);
            C(tic,i)=cohesion(3);
            n=n+1;
            if tally==0
                tic=tic+1;
            end
        end
    end
end

% round thicknesses to prevent zero layer thickness error
thickness=round(10*thickness)/10;
    
% Hanson and Simon, 2001 conversion for Kb using taucrit
kb=0.2*C.^-0.5;

% Hack for multiple grain partitions test. Each fault layer gets a specific
% fraction value
gfrac1=99.99999;
gfrac2=90;
gfrac3=10;
gfrac4=00.00001;

% Proceeding section edited from Greg's code for making a layer file.
%--now, make the lay file--%
% Create layer file and open it for writing
filenm= [ basenm '.lay' ];
lfid = fopen(filenm,'w');
if lfid<1
    error('Unable to create layer file.\n');
end

fprintf('now writing lay file.\n');

% Write header information
fprintf(lfid,' %.2f\n',ts);
fprintf(lfid,'%d\n',length(xc));

% Loop over nodes to write information
for j=1:length(xc)
    nlayers(j)=sum(thickness(:,j)>0);
    layindex=find(thickness(:,j));
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',nlayers(j));
    
    % For each layer at this node, write information
    for i=1:nlayers(j)
        % Basic properties: creation time, recent activity time, exposure
        % time, thickness, erodibility, and regolith/bedrock flag
        fprintf(lfid,'%.2f %.2f %.2f\n',0,0,0);
        %         if C(layindex(i),j)==cohesion(4)
        fprintf(lfid,'%.2f %f %.0f\n',thickness(layindex(i),j),kb(layindex(i),j),0); %end
        %         if C(layindex(i),j)==cohesion(3) || C(layindex(i),j)==cohesion(2)
        %             fprintf(lfid,'%.2f %f %.0f\n',thickness(layindex(i),j),kb(layindex(i),j),1); end
        if numg==1
            fprintf(lfid,'%.2f\n',thickness(layindex(i),j));
        end
%         if numg==2
%             if C(layindex(i),j)==cohesion(4);
%                 fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac1),(thickness(layindex(i),j)*gfrac2)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
%             if C(layindex(i),j)==cohesion(3) || C(layindex(i),j)==cohesion(2) || C(layindex(i),j)==cohesion(1)
%                 fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac3),(thickness(layindex(i),j)*gfrac4)); end
%         end
        if numg==2
            if C(layindex(i),j)==cohesion(4);
                fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac1),(thickness(layindex(i),j)*gfrac4)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
            if C(layindex(i),j)==cohesion(3);
                fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac2),(thickness(layindex(i),j)*gfrac3)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
            if C(layindex(i),j)==cohesion(2);
                fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac3),(thickness(layindex(i),j)*gfrac2)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
            if C(layindex(i),j)==cohesion(1);
                fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac4),(thickness(layindex(i),j)*gfrac1)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
            if C(layindex(i),j)==cohesion(3) || C(layindex(i),j)==cohesion(2) || C(layindex(i),j)==cohesion(1)
                fprintf(lfid,'%.2f %.2f\n',(thickness(layindex(i),j)*gfrac3),(thickness(layindex(i),j)*gfrac4)); end
        end

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

