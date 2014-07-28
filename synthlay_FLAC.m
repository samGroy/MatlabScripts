function product = synthlay_FLAC( basenm,b,cc,numg, gfrac1, gfrac2, gfrac3, gfrac4)
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

% SGR 10/2013 THOUGHT
% Current coupled code employs kb map, not layfiles. Though it might be
% slightly slower, if we switch back to use of layer files we could control
% grain size, and therefore have coupled transport limited models. This
% could make lots of people happy. "Layers" could be vertically continuous.
    
cc=cc(b==0);
% Hanson and Simon, 2001 conversion for Kb using taucrit
kb=0.2*cc.^-0.5;
ccmean=1E7;
% Thickness, 10km should do the trick
thickness=10000;
% Hack for multiple grain partitions test. Each fault layer gets a specific
% fraction value
gfrac1=.99999;
% gfrac2=80;
% gfrac3=20;
gfrac4=.00001;

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
fprintf(lfid,' %.0f\n',0);
fprintf(lfid,'%d\n',length(cc));

% Loop over nodes to write information
for j=1:length(cc)
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',1);
    
    % For each layer at this node, write information
    % Basic properties: creation time, recent activity time, exposure
    % time, thickness, erodibility, and regolith/bedrock flag
    fprintf(lfid,'%.2f %.2f %.2f\n',0,0,0);
    %         if C(layindex(i),j)==cohesion(4)
    fprintf(lfid,'%.2f %f %.0f\n',thickness,kb(j),0); %end
    %         if C(layindex(i),j)==cohesion(3) || C(layindex(i),j)==cohesion(2)
    %             fprintf(lfid,'%.2f %f %.0f\n',thickness(layindex(i),j),kb(layindex(i),j),1); end
    if numg==1
        fprintf(lfid,'%.2f\n',thickness(j));
    end
    if numg==2
        if cc(j)>=ccmean;
            fprintf(lfid,'%.2f %.2f\n',(thickness*gfrac4),(thickness*gfrac1)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
        if cc(j)<=ccmean;
            fprintf(lfid,'%.2f %.2f\n',(thickness*gfrac1),(thickness*gfrac4)); end %take gfrac 1 and 3 as the smaller grain size: proportion 1
    end
end
fclose(lfid);
fprintf('done.\n');

