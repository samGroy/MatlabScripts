% TEST 4: tighter loop architecture, see if there is speedup if everything is saved to ramdisk.
% Result: ~5% speedup before regrid, unknown after...~1200 ka in ~9 hours
% 9/20/2013 essentially halves sim time in total (beyond 1 Ma)
%runflacchild: couples FLAC and CHILD
%
% 12/2012 SGR In prep for the 2013 geomorph grant proposal, cohesion
% transfer from FLAC is now online, uses kb map, so vertically continuous.
% CHILD can handle fine mesh, but FLAC can't, so it is always the first to
% give up.
% 6/13 SGR lateral displacement is now transferred from FLAC to CHILD. This
% is done by using a points file input. The disadvantages are 1) CHILD will
% not keep track of cumulate time, 2) we can't define cohesion in 3D
% (although we haven't done so anyway). Advantages are 1) lateral movement
% is possible, 2) heterogeneous grid resolution is also possible, so we can
% refine the mesh around interesting locations.
% INITIALIZE

% 1) Assign dirs, define choice CHILD parameters
%   LOOP
% 2) Run FLAC,if first time, find surface dimensions, assign boundary codes
% 3) convert flac dat for child readability (including randomization of xyz)
% 4) Run CHILD for dt
% 5) convert child elevation data for flac input

% Q: is it faster to write new small files every run (flacdat) or append data to files that already exist?
% Q: do I want to append all child data into their respective files (tri,
% nodes, z, etc) or make a superfile that contains all data needed? Or keep
% going with the flacdat files, they make for good videos...

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Switches for implementing cohesion and saving a series of the flacdat
%files
opt_ramdisk=1; %switch for ramdisk, really worthwhile as it cuts massive time.
opt_dynamic_erody = 2;  % Option to use flac cohesion field for child erodibility
opt_numbered_flac_files = 1; % keep record of previous flac outputs
opt_continue = 0; %option to continue a run, number here equals the iteration, itr, you want to continue from
opt_boundaries = 2; %how many open boundaries
opt_child_output = 1; %output child data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if opt_continue
    start=opt_continue;
else
    %%%%% Initialization %%%%%
    fprintf('Flacchild simulation in progress\n...\n');
    % User-adjustable parameters and names here:
    % General settings:
    mysys='dos';  % Phaedra: change this to dos
    filesys=[''];
    
    % Time parameters
    child_time = 0;
    Dt = 1000;  % Duration for each sub-run
    runtime = 2000000;  % Total desired run time in years, FLAC will usually fail before long time durations complete
    nt = runtime/Dt; %number of time steps
    next_child_output = Dt;
    start=1;
    
    %Directories, file names
    % dos('imdisk -a -s 3000M -m X: -p "/fs:fat32 /q /y"'); %Make sure ramdisk
    % is mounted!
    
    loaddir = 'C:\child_n\ChildExercises\FLACCHILD\';
    if opt_ramdisk==1
        rundir = 'X:\'; % where CHILD and this script live. Must run imdisk before this script!!
    else
        dos(['mkdir ' loaddir '\newrun\']);
        rundir=([loaddir '\newrun\']);
    end
    % savdir = '\child_n\ChildExercises\gregsFlacChild\newrun\'; %where output files are saved
    % Initialize names and places if the directory does not exist
    % dos(['mkdir ' savdir]); % create the save folder
    cd(loaddir)
    
    
    % FLAC stuff
    flacinputfilenm = 'flacdat';% (base) name of flac output file
    flacinputfullname = [rundir flacinputfilenm];
    flacsetupfile = 'setup_flac_child.txt'; % Name of initial flac script
    flacrunfile = 'run_flac_child.txt'; % Name of flac run script
    flacfilepath = '"C:\Program Files\Itasca\Flac3d400\exe64\flac3d400_gui_64"'; % where flac lives
    if opt_numbered_flac_files
        flacextension = '.txt';
        flacfilenumber = 1;
        current_flac_file = [flacinputfullname num2str(flacfilenumber) flacextension];
    else
        flacfilenm = [rundir flacinputfilenm];
        current_flac_file = flacinputfullname;
    end
    strainswitch=0; %for opt_dynamic_erody=2 when plastic failure initiates
    stweak=0.03;
    old_pks=[];
    
    % Child stuff:
    childexe = [rundir 'child'];
    childptsfilenm = 'childgridconfig2.pts'; % Points interpolated from FLAC used to create CHILD surface
    kbfile = 'kb.txt'; % Interpolated erodibility for CHILD
    childinputfiletemplate = 'cftemplate.in'; % input file template
    childrunbasenm = 'testflacchild'; %name of CHILD model run
    childdat = 'childdat';
    childinputfile = [childrunbasenm '.in'];
    r=20; % magnitude of mesh refinement from FLAC to CHILD
    % Get a copy of child's input file so we can modify it
    
    ciftemplate = creadinfile([loaddir childinputfiletemplate]); %read from loaddir...
    cif = ciftemplate;
    % Edit Input file for initial CHILD run.
    cif = cchangeparam(cif,'RUNTIME',Dt);
    cif = cchangeparam(cif,'OPINTRVL',Dt);
    cif = cchangeparam(cif, 'ST_PMEAN',.0001); %SGR kept at 0.01 for buildup model
    cif = cchangeparam(cif, 'ST_STDUR',10);
    cif = cchangeparam(cif, 'ST_ISTDUR',0);
    cif = cchangeparam(cif, 'OPTREADINPUT',12);
    cif = cchangeparam(cif,'OPTLAYEROUTPUT',0);
    if opt_dynamic_erody %== 1
        cif = cchangeparam(cif,'OPT_SET_ERODY_FROM_FILE',1);
        cif = cchangeparam(cif,'ERODYFILE_NAME',kbfile);
        % else
        %     cif = cchangeparam(cif,'OPT_SET_ERODY_FROM_FILE',0);
        %     cif = cchangeparam(cif,'KB',0);
    end
    cwriteinfile(childinputfile,cif); %saves edited template as a new in file SGR
    
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Move stuff to ramdisk X 9/16/2013
    [s,w]=dos(['del ' rundir childdat '.*']); % Clear ramdisk of child files before adding files, make sure things are saved elsewhere if needed!
    % [s,w]=dos(['copy ' loaddir flacsetupfile ' ' rundir]);
    [s,w]=dos(['copy ' loaddir flacrunfile ' ' rundir]);
    [s,w]=dos(['copy ' loaddir 'REGRID1.txt ' rundir]);
    % [s,w]=dos(['copy ' loaddir 'flac_ini.f3sav ' rundir]); %omit if running from scratch or save file already positioned
    [s,w]=dos(['copy ' loaddir 'child.exe ' rundir]);
    [s,w]=dos(['copy ' loaddir 'testflacchild.in ' rundir]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%% Cycle %%%%%%
cd(rundir)
for itr=start:nt
    fprintf('time: %.f ka\n',itr);
    % LEAD CYCLE WITH FLAC SETUP RUN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Run the FLAC setup file, only run once, rerun if the setup file is changed
    if itr==1
        [s,wf]=dos(['start /min /wait "" ' flacfilepath ' call ' flacsetupfile ' quit']);% initiate flac run%%%%%%%%%%%%%%%%%%%%%%%%
    else
        [s,wf]=dos(['start /min /wait "" ' flacfilepath ' call ' flacrunfile ' quit']);% initiate flac run%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if s>0
        fprintf('Error message from FLAC:\n\n%s',wf);
        error('Bailing out of runflacchild.m');
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % Rename the next flac data file
    if opt_numbered_flac_files
        flacfilenumber = itr;
        current_flac_file = [rundir flacinputfilenm num2str(flacfilenumber) '.txt'];
    end
    
    % Read the FLAC data
    [idf,xf,yf,zf,vf,cf,ef,pf]=readflac2(current_flac_file);
    
    
    % Find mesh boundaries, assign boundary code, apply slight randomization to child nodes
    % first cycle only
    if itr==1
        Lx = max(xf);    % Grid length in m
        Ly = max(yf);    % Grid width in m
        dxf=max(diff(xf)); dyf=max(diff(yf)); % x and y grid spacing
        %         dxc=dxf/r; dyc=dyf/r; % Child resolution is Xx the flac resolution
        xc_res=1500; yc_res=xc_res; %hack to make nonuniform grid uniform and fine. MUST BE FACTOR OF Lx or you lose the west or north boundary!!!!
        [x,y]=meshgrid(min(xf):xc_res:Lx,min(yf):yc_res:Ly); %get coordinate values
        x=x'; y=y'; %It is necessary to rotate for the next step...
        xc=reshape(x,1,[]); % Flatten x y matrices to vector format
        yc=reshape(y,1,[]);
        F=TriScatteredInterp(xf',yf',zf'); %interp to get elevation data
        zc=F(xc,yc); % z as function of x,y thanks to flac data
        
        %     Boundary code assumes mesh is rectangular
        b = zeros(size(xc));   % boundary codes: initially all zeros
        if opt_boundaries == 4
            b( xc<xc_res/2 ) = 2; % all nodes with x coordinate less than half the node spacing, flac grid is uniform
            b( xc>(Lx-xc_res/2) ) = 2; %other side, etc... all open boundaries
            b( yc<yc_res/2 ) = 2;
            b( yc>(Ly-yc_res/2) ) = 2;
        else
            b( xc<xc_res/2 ) = 1; % closed east and west boundaries
            b( xc>(Lx-xc_res/2) ) = 1; 
            b( yc<yc_res/2 ) = 2; % open north and south boundaries
            b( yc>(Ly-yc_res/2) ) = 2;
        end
        child_interior_nodes = find( b==0 );
        
        %
        % Create child's node coordinates by adding small random offsets to each
        % node
        % >varres: need to find out number of child nodes versus FLAC nodes
        random_offset_x = zeros(1,length(xc));  % slightly offset flac coords, keep track of the offset so you can keep track between flac and child nodes later
        random_offset_y = zeros(1,length(xc));  % slightly offset flac coords
        random_offset_z = zeros(1,length(xc));  % slightly offset flac coords, z not really necessary unless you want random surface for more natural discharge channeling SGR 12/19/12
        random_offset_x(child_interior_nodes) = (xc_res/4)*rand(1,length(child_interior_nodes));  % slightly offset flac coords
        random_offset_y(child_interior_nodes) = (yc_res/4)*rand(1,length(child_interior_nodes));  % slightly offset flac coords
        random_offset_z(child_interior_nodes) = rand(1,length(child_interior_nodes));  % slightly offset flac coords
        xc = xc+random_offset_x;  % slightly offset flac coords
        yc = yc+random_offset_y;
        zc = zc+random_offset_z;
    else
        %         if opt_dynamic_erody == 2
        %             [xc,yc,s,b,zc]=creadxyzb2(childrunbasenm,2);
        %         end
        % Find flac's offset in x, y, and z, and interpolate to child
        dxf = xf - xf_prev;
        dyf = yf - yf_prev;
        dzf = zf - zf_prev;
        %     dxc=griddata(xf+tiny_offset_x,yf+tiny_offset_y,dxf,xc,yc); Uncomment
        %     if CHILD chokes
        %     ddxc=griddata(xf+tiny_offset_x,yf+tiny_offset_y,dxf,xc,yc,'nearest');
        %     dxc(isnan(dxc))=ddxc(isnan(dxc));  % in case corner points come out as NaN
        %     dyc=griddata(xf+tiny_offset_x,yf+tiny_offset_y,dyf,xc,yc);
        %     ddyc=griddata(xf+tiny_offset_x,yf+tiny_offset_y,dyf,xc,yc,'nearest');
        dxc=griddata(xf,yf,dxf,xc,yc);
        ddxc=griddata(xf,yf,dxf,xc,yc,'nearest');
        dxc(isnan(dxc))=ddxc(isnan(dxc));  % in case corner points come out as NaN
        dyc=griddata(xf,yf,dyf,xc,yc);
        ddyc=griddata(xf,yf,dyf,xc,yc,'nearest');
        dyc(isnan(dyc))=ddyc(isnan(dyc));  % in case corner points come out as NaN
        Fflac = TriScatteredInterp(xf',yf',dzf');
        Fflac_nearest  = TriScatteredInterp(xf',yf',dzf','nearest');
        dzc = Fflac(xc,yc);
        dzcn = Fflac_nearest(xc,yc);
        dzc(isnan(dzc))=dzcn(isnan(dzc));
        
        xc = xc + dxc;
        yc = yc + dyc;
        zc = zc + dzc;
        
        for i=1:length(b) %this sets the boundary elevation to be equal to its closest interior neighbor SGR 1/10/13
            if b(i)>0
                dx2=xc-xc(i); dx2=dx2.^2;
                dy2=yc-yc(i); dy2=dy2.^2;
                dist=sqrt(dx2+dy2);
                dist(b>0)=100000000;
                [val,ind]=min(dist);
                zc(i)=zc(ind);
            end
        end
    end
    
    %     MAP COORDINATES AND ERODIBILITY
    if opt_dynamic_erody == 1 || max(pf)<=stweak
        
        %FIRST: surface moved, make a points file for CHILD
        fid=fopen([rundir childptsfilenm],'w');
        if fid<=0, error('Unable to create file for CHILD points');end
        fprintf(fid,' %d\n',length(xc));  % print # of points in file
        for i=1:length(xc) % algo writes xyzb column data in file
            fprintf(fid,'%.4f %.4f %.3f %d\n',xc(i),yc(i),zc(i),b(i) );
        end
        fclose(fid);
        
        %LAST: make kb map for CHILD directly from FLAC data
        cc=griddata(xf,yf,cf,xc,yc);
        ccn=griddata(xf,yf,cf,xc,yc,'nearest');
        cc(isnan(cc))=ccn(isnan(cc)); %SGR remove nans 1/2/13
        %         cc(cc<=2.9E7)=1E5; %SGR added this weak/strong division 11-13-2013
        cerody=(0.2*cc.^-.5); %SGR 1/2/13 Hanson and Simon '01 equation
        %         cerody(xc<=Lx/2)=cerody(xc<=Lx/2)*10; %just added this for side-by-side comparison 5/1/14
        fid=fopen([rundir kbfile],'w');
        if fid<=0, error('Unable to create file for CHILD points');end
        for i=1:length(cerody) % algo writes xyzb column data in file
            fprintf(fid,'%f\n',cerody(i)); %consider changing to weak/strong division
        end
        fclose(fid);

        % synch FLAC weak zones to high resolution CHILD damage zones
    elseif opt_dynamic_erody == 2 && max(pf)>=stweak 
        if strainswitch==0
            fprintf('Permanent crustal strain initiated!\n');
            strainswitch=1;
            faultfile=zeros(1,7);
            c1=5e5; c2=1e7; c3=3e7-.1; c4=3e7;
            faultfile(1,1)=c1; faultfile(1,2)=c2; faultfile(1,3)=c3; faultfile(1,4)=c4;
            thick=yc_res;
            strike=0; dip=30;
            
            ciftemplate = creadinfile([rundir childinputfile]); %read from loaddir...
            cif = ciftemplate;
            cif = cchangeparam(cif,'OPT_SET_ERODY_FROM_FILE',0);
            cif = cchangeparam(cif,'OPTREADINPUT',1);
            cif = cchangeparam(cif,'INPUTDATAFILE',childrunbasenm);
            cif = cchangeparam(cif,'INPUTTIME',Dt);
            cif = cchangeparam(cif,'OPTREADLAYER',1);
            cif = cchangeparam(cif,'OPTLAYEROUTPUT',1);
            cwriteinfile(childinputfile,cif); %saves edited template as a new in file SGR
        end
        
        %FIRST: surface changed, create new node and z files
        fid=fopen([rundir childrunbasenm '.nodes'],'w');
        if fid<=0, error('Unable to create file for CHILD nodes');end
        fprintf(fid,' %d\n',Dt);  % print the time
        fprintf(fid,'%d\n',length(xc));  % print # of points in file
        for i=1:length(xc) % algo writes xyzb column data in file
            fprintf(fid,'%.4f %.4f %d %d\n',xc(i),yc(i),sp(i),b(i) );
        end
        fclose(fid);
        fid=fopen([rundir childrunbasenm '.z'],'w');
        if fid<=0, error('Unable to create file for CHILD elev');end
        fprintf(fid,' %d\n',Dt);  % print the time
        fprintf(fid,'%d\n',length(zc));  % print # of points in file
        for i=1:length(zc) % algo writes xyzb column data in file
            fprintf(fid,'%.4f \n',zc(i));
        end
        fclose(fid);
        
        %LAST: plastic strain initiated, find localizations and plant
        %damage zones on them
        yline=yf(xf<100 & cf<c4); pline=pf(xf<100 & cf<c4); zline=zf(xf<100 & cf<c4);
        [pks,locs]=findpeaks(pline);
        fprintf('%d faults located\n',length(pks));
        if length(pks)==length(old_pks) % IF THE SAME STRAIN PEAKS EXIST, WITH NO NEW ADDITIONS
            system(['rm ' childrunbasenm '.lay']);
            system(['ren ' childrunbasenm '.lay1 ' childrunbasenm '.lay']); %recycle the old layer file
        else %IF NEW STRAIN PEAKS APPEAR, THEY NEED TO BE ADDED TO THE FAULTFILE
            old_pks=pks;
            ycoords=yline(locs); zcoords=zline(locs);
            for i=1:length(pks)
                if ycoords(i)>Ly/2
                    faultfile(i+1,1)=dip;
                else
                    faultfile(i+1,1)=-dip;
                end
                faultfile(i+1,2)=strike;
                faultfile(i+1,3)=thick;
                faultfile(i+1,4)=1;
                faultfile(i+1,5)=Lx/2;
                faultfile(i+1,6)=ycoords(i);
                faultfile(i+1,7)=zcoords(i);
            end
            synthlay_grad7( childrunbasenm, faultfile, Dt, 1, 1, 1, 1, 1, 1) % Should automatically spit out a .lay file for use
        end
    end
    
    
    
    % Run Child
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    runcommand = [childexe ' ' rundir childrunbasenm '.in --no-check'];
    % Now execute the run
    if strcmp(mysys,'unix')
        [s,w]=unix(runcommand);
    elseif strcmp(mysys,'dos')
        [s,w]=dos(runcommand);
    end
    if s>0
        fprintf('Error message from CHILD:\n\n%s',w);
        error('Bailing out of runflacchild.m');
    end
    %\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
    
    %save the CHILD surface in a struct that can be viewed after running
    %load tri, x, y, z data
    %         h=ctrisurf(childrunbasenm,1);
    %         m(itr)=getframe;
    %
    % Advance the "child clock"
    child_time = child_time + Dt;
    
    % Get the new child x,y,z coords
    xc_prev = xc;  % remembering the old coords isn't strictly necessary
    yc_prev = yc;  % ...but it helps w/ debugging
    zc_prev = zc;
    
    [xc,yc,sp,b,zc]=creadxyzb2(childrunbasenm,2);% Get child coords and node data
    
    %     Interpolate child's elevation field to flac nodes
    F=TriScatteredInterp(xc,yc,zc);
    Fn=TriScatteredInterp(xc,yc,zc,'nearest');
    zf = F(xf,yf);
    zfn=Fn(xf,yf);
    zf(isnan(zf))=zfn(isnan(zf));  % in case corner points come out as NaN
    
    % Write interpolated elevation data to a file that flac can understand
    flacoutfid = fopen('topography4flac.txt','wt'); % Open output file
    for i=1:length(zf)
        fprintf(flacoutfid,'%e\n',zf(i));
    end
    fclose(flacoutfid);
    
    %     % Remember flac's old coordinates before rerunning
    xf_prev = xf;
    yf_prev = yf;
    zf_prev = zf;
    
    %if you asked for it, save the child data to a text file
    if opt_child_output
        if strainswitch==1
            k=readlaykb3(childrunbasenm,2,1,length(xc));
            k=k';
        else
            if itr==1
            k=log10(cerody);
            else
            k=log10(cerody');
            end
        end
        xyzk=[xc'; yc'; zc'; k];
        save('childdat.txt', 'xyzk', '-ASCII','-append');
    end
    %     Append new CHILD data to previous data
%     system(['type ' childrunbasenm '.nodes ' '>>' childdat '.nodes']); %appends new data to old list
%     system(['type ' childrunbasenm '.z ' '>>' childdat '.z']); %appends new data to old list
%     system(['type ' childrunbasenm '.q ' '>>' childdat '.q']); %appends new data to old list
%     system(['type ' childrunbasenm '.area ' '>>' childdat '.area']); %appends new data to old list
%     system(['type ' childrunbasenm '.net ' '>>' childdat '.net']); %appends new data to old list
%     system(['type ' childrunbasenm '.slp ' '>>' childdat '.slp']); %appends new data to old list
%     system(['type ' childrunbasenm '.tau ' '>>' childdat '.tau']); %appends new data to old list
%     system(['type ' childrunbasenm '.varea ' '>>' childdat '.varea']); %appends new data to old list
%     system(['type ' childrunbasenm '.tri ' '>>' childdat '.tri']); %appends new data to old list
    
    % If it's time, save child files to a special folder
    %     if child_time>=next_child_output || itr==nt
    %         myfoldername = [childrunbasenm num2str(itr)];
    %         if strcmp(mysys,'unix')
    %             unix(['mkdir ' myfoldername]);
    %             mycommand = ['cp ' childrunbasenm '.* ' myfoldername];
    %             [s,w]=unix( mycommand );
    %         elseif strcmp(mysys,'dos')
    %             dos(['mkdir ' myfoldername]);
    %             [s,w]=dos(['copy ' childrunbasenm '.* ' myfoldername]);
    %         else
    %             error('You need to set mysys to either unix or dos');
    %         end
    %         next_child_output = next_child_output + 10*Dt;
    %     end
    fclose all;
end

