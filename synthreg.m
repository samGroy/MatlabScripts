function  synthreg( basenm,step)
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
% Added input functions for gran size info.

regthick=5;
nlayers=2;
kb=([1E-1 1E-5]);

cname = basenm;
xyzb = creadxyzb(cname,step);
nodes=length(xyzb);
for i=1:nodes
    if xyzb(i,4)==0
        xc(i) = xyzb(i,1);
        yc(i) = xyzb(i,2);
        zc(i) = xyzb(i,3);
        zci(i)=zc(i);
        bc(i) = xyzb(i,4);
    end
end
intnode=length(xc);
%--now, make the lay file--%
% Create layer file and open it for writing (first make sure it doesn't
% already exist)

filenm= [ basenm '.lay' ];

    %     filenm= [ basenm sprintf('.lay%.0f',step) ];
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


% Write header information
fprintf(lfid,' %.2f\n',0);
fprintf(lfid,'%d\n',length(intnode));

% Loop over nodes to write information
for j=1:length(xc)
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',nlayers);
    
    % For each layer at this node, write information
    for i=1:nlayers
        % Basic properties: creation time, recent activity time, exposure
        % time, thickness, erodibility, and regolith/bedrock flag
        fprintf(lfid,'%.2f %.2f %.2f\n',0,0,0);
        %         if C(layindex(i),j)==c1
        if i==1
            fprintf(lfid,'%.2f %f %.0f\n',regthick,kb(1),0); %end
            fprintf(lfid,'%.2f\n',regthick);
        end
        if i==2
            fprintf(lfid,'%.2f %f %.0f\n',regthick*10000,kb(2),0); %end
            fprintf(lfid,'%.2f\n',regthick*10000);
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



