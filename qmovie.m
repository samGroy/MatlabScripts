function [m,data,erate]=qmovie(basenm,nframes,ac,az,numg,ax,ay)
% CMOVIE: animation of CHILD surface (or any other field) through time.
%    Usage: m = cmovie( basenm, nframes, ax, ay, az, ac )
%
% Parameters:
%   basenm - base file name of CHILD run
%   nframes - number of frames
%   ax (optional) - extent of x-axis
%   ay (optional) - extent of y-axis
%   az (optional) - extent of z-axis
%   ac (optional) - parameter code: 1=elevation, 2=discharge, 3=channel
%   width, 4=slope, 5=shear stress, 6=uplift, 7=total erosion per node,
%   8=tortuosity at node spacing resolution,9=finest grain size fraction,
%   10=kb, 11=change in slope over timestep, 12=calculate erosion rate over
%   timesteps for a determined spatial range, 13=precip, 14=drainage area,
%   15=sed flux, 16=sed influx from fluvial, 17=sed influx from diffusion,
%   18=regolith thickness, 19=bedrock/regolith flagmap (% coverage time), 20=total layer
%   exposure time, 21=layer recent activity time, 22=layer creation time,
%   23=number of layers, 24=cumulative bedrock exposure time, 25=dzdt
%
%

% Open files
g=1;
time=0;
average=0;
flacchild=0;
filesys='';
filenm= [filesys basenm '.nodes' ];
nfid=fopen(filenm,'r');
if nfid<=0, error('Unable to open node file'),end
filenm= [filesys basenm '.tri' ];
tfid=fopen(filenm,'r');
if tfid<=0, error('Unable to open triangle file'),end
filenm= [filesys basenm '.z' ];
zfid=fopen(filenm,'r');
if zfid<=0, error('Unable to open elevation file'),end
if ac==2
    fprintf('reading channel discharge data \n');
    filenm=[filesys basenm '.q'];
    qfid=fopen(filenm,'r');
end
if ac==3
    fprintf('reading channel width data \n');
    filenm=[filesys basenm '.chanwid'];
    wfid=fopen(filenm,'r');
    if wfid<=0
        fprintf('empirical width equation used \n');
        filenm=[filesys basenm '.q'];
        qfid=fopen(filenm,'r');
    end
    if wfid>0
        fprintf('Finnegan width equation used \n');
    end
end
if ac==4
    fprintf('reading slope data \n');
    filenm=[filesys basenm '.slp'];
    sfid=fopen(filenm,'r');
end
if ac==5
    fprintf('reading basal shear stress data \n');
    filenm=[filesys basenm '.tau'];
    taufid=fopen(filenm,'r');
end
if ac==6 || ac==12
    fprintf('reading uplift data \n');
    filenm=[filesys basenm '.up'];
    ufid=fopen(filenm,'r');
end
if ac==7
    fprintf('calculating erosion \n');
    filenm= [filesys basenm '.z' ];
    zfid=fopen(filenm,'r');
    filenm=[filesys basenm '.up'];
    ufid=fopen(filenm,'r');
%     filenm= [filesys basenm '.dzdt' ];
%     zfid=fopen(filenm,'r');
end
if ac==9
    fprintf('reading grain size fraction data \n');
    filenm= [filesys basenm '.tx' ];
    txfid=fopen(filenm,'r');
end
% fhandle=figure;
if ac==11
    fprintf('reading slope data \n');
    filenm=[filesys basenm '.slp'];
    sfid=fopen(filenm,'r');
end
if ac==13
    fprintf('reading flux and area data \n');
    filenm=[filesys basenm '.q'];
    qfid=fopen(filenm,'r');
    filenm=[filesys basenm '.area'];
    afid=fopen(filenm,'r');
end
if ac==14
    fprintf('reading area data \n');
    filenm=[filesys basenm '.area'];
    afid=fopen(filenm,'r');
end
if ac==15
    fprintf('reading sed flux data \n');
    filenm=[filesys basenm '.qs'];
    qsfid=fopen(filenm,'r');
end
if ac==16
    fprintf('reading sed flux data \n');
    filenm=[filesys basenm '.qsin'];
    qsfid=fopen(filenm,'r');
end
if ac==17
    fprintf('reading sed flux data \n');
    filenm=[filesys basenm '.qsdin'];
    qsfid=fopen(filenm,'r');
end
if ac==25
    fprintf('reading dzdt data \n');
    filenm=[filesys basenm '.dzdt'];
    dzdtfid=fopen(filenm,'r');
end
for i=1:nframes
    
    % Read stuff
    tm = fscanf(nfid,'%f',1);
    fprintf('CTRISURF: Reading time %f\n',tm);
    tm = fscanf(tfid,'%f',1);
    nn = fscanf(nfid,'%d',1);
    nt = fscanf(tfid,'%d',1);
    t=fscanf(tfid,'%f',[9,nt]);
    n=fscanf(nfid,'%f',[4,nn]);
    tm = fscanf(zfid,'%f',1);
    nn = fscanf(zfid,'%d',1);
    z=fscanf(zfid,'%f',[1,nn]);
    data=zeros(3,length(n));
    data(1:2,:)=n(1:2,:);
    data=data';
    
    if ac==1 || ac==12, c=z; text='elevation'; end
    if ac==2
        tm = fscanf(qfid,'%f',1);
        nn = fscanf(qfid,'%d',1);
        c=fscanf(qfid,'%f',[1,nn]);
        c(z<=0)=0;
        if average==1
            if i==1, flagmatrix=zeros(nn,1); end
            flagmatrix(:,i)=c;
            c=mean(flagmatrix,2);
        end

        text='discharge';
    end
    if ac==3
        if wfid>0
            tm = fscanf(wfid,'%f',1);
            nn = fscanf(wfid,'%d',1);
            c=fscanf(wfid,'%f',[1,nn]);
        end
        if wfid<=0
            tm = fscanf(qfid,'%f',1);
            nn = fscanf(qfid,'%d',1);
            q=fscanf(qfid,'%f',[1,nn]);
            secperyear=31556926;
            c=10*(q/secperyear).^.5;
        end
        text='channel width';
    end
    if ac==4
        tm = fscanf(sfid,'%f',1);
        nn = fscanf(sfid,'%d',1);
        c=fscanf(sfid,'%f',[1,nn]);
        text='slope';
    end
    if ac==5
        tm = fscanf(taufid,'%f',1);
        nn = fscanf(taufid,'%d',1);
        c=fscanf(taufid,'%f',[1,nn]);
        text='shear stress';
%         if mean(diff(c))==0, continue, end %  commin for coupled runs
%         only
    end
    if ac==6 || ac==12
        tm = fscanf(ufid,'%f',1);
        nn = fscanf(ufid,'%d',1);
        if ac==12
            c2=rot90(fscanf(ufid,'%f',[1,nn]));
        else
            c=fscanf(ufid,'%f',[1,nn]);
            text='uplift rate';
        end
    end
    if ac==7
        % zed=z;
        % tm = fscanf(ufid,'%f',1);
        % nn = fscanf(ufid,'%d',1);
        % up=fscanf(ufid,'%f',[1,nn]);
        % c=tm.*up-zed;
        if i==1
            c=z;
            text='erosion';
            tmold=tm;
            zold=z;
%         elseif mod(i,2)==0 % for flacchild runs
% %         elseif i>1
%             c=(zold-z)/(tm-tmold);
% %             c=(zold-z)/1000; %for flacchild runs 10-9-2013
%             %       text='total erosion';
            text='erosion';
        else
            c=(z-zold)./(tm-tmold);
            zold=z;
            tmold=tm;
%             tmold=tm;
%             zold=z;
%         elseif mod(i,2)==1, zold=z; continue
        end
    end
    if ac==8
        c=tortuositymap(basenm,i);
        text='tortuosity at node resolution';
    end
    if ac==9
        tm = fscanf(txfid,'%f',1);
        nn = fscanf(txfid,'%d',1);
        c=fscanf(txfid,'%f',[1,nn]);
        if average==1
            if i==1, flagmatrix=zeros(nn,1); end
            sedflag=readlaykb2(basenm,i,numg,nn,19);
            c(sedflag==0)=nan;
            flagmatrix(:,i)=c;
            c=nanmean(flagmatrix,2);
            c(isnan(c))=100;
%             scatter3(n(1,:),n(2,:),z,10,c);
%             for ii=1:nn
%                 c(ii)=mean(flagmatrix(ii,(flagmatrix(ii,:)>=0)),2);
%             end
        end
        text='smallest grain size fraction';
    end
    if ac==10 || ac==10.1 || ac==18 || ac==20 || ac==21 || ac==22 || ac==23 || ac==24  || ac==19
        if average==0
            if exist('numg','var')==0, numg=1; end
            c=readlaykb2(basenm,i,numg,nn,ac);
            text=(['layer data ' num2str(ac)]);
        elseif average==1
            if i==1, flagmatrix=zeros(nn,1); end
            [c,flagmatrix]=readlaymean(basenm,i,ac,numg,nn,flagmatrix);
            text=(['layer data ' num2str(ac)]);
        end
    end
    if ac==11
        if i==1
            tm = fscanf(sfid,'%f',1);
            nn = fscanf(sfid,'%d',1);
            c=fscanf(sfid,'%f',[1,nn]);
            text='slope';
            tmold=tm; sold=c;
        else
            tm = fscanf(sfid,'%f',1);
            nn = fscanf(sfid,'%d',1);
            s=fscanf(sfid,'%f',[1,nn]);
            c=(s-sold)/(tm-tmold);
            text='slope';
            tmold=tm; sold=s;
        end
    end
    if ac==13 %works kind of but doesn't work for streams that pass from large rainfall area to small rainfall area.
        tm = fscanf(qfid,'%f',1);
        nn = fscanf(qfid,'%d',1);
        qarray=fscanf(qfid,'%f',[1,nn]);
        tm = fscanf(afid,'%f',1);
        nnarea = fscanf(afid,'%d',1);
        areaarray=fscanf(afid,'%f',[1,nnarea]);
        c=qarray(1:nnarea)./areaarray(1:nnarea);
        c(nnarea+1:nn)=0;
        text='precipitation rate';
    end
    if ac==14
        tm = fscanf(afid,'%f',1);
        nnarea = fscanf(afid,'%d',1);
        c=fscanf(afid,'%f',[1,nnarea]);
        c(nnarea+1:nn)=0;
        text='drainage area';
    end
    if ac==15 || ac==16 || ac==17
        tm = fscanf(qsfid,'%f',1);
        nn = fscanf(qsfid,'%d',1);
        c=fscanf(qsfid,'%f',[1,nn]);
        c(nn+1:nn)=0;
        if average==1
            if i==1, flagmatrix=zeros(nn,1); end
            flagmatrix(:,i)=c;
            c=mean(flagmatrix,2);
        end
        text='sed flux';
    end
    if ac==25
        tm = fscanf(dzdtfid,'%f',1);
        nn = fscanf(dzdtfid,'%d',1);
        c=fscanf(dzdtfid,'%f',[1,nn]);
        text='dzdt';
    end
    
    
    % Extract coordinates
    tri = [ rot90(t(1,:),3) rot90(t(2,:),3) rot90(t(3,:),3)]+1;
    x = n(1,:);
    y = n(2,:);
    b = n(4,:);
    
    % Remove edge effects
    z = cinterpclosededges( x, y, z, b );  % this removes edge effects for plotting
    
    % Plot
    h=trisurf(tri,x,y,z,c);%,'EdgeAlpha',0.1); %edgealpha crashes matlab
    set(h,'FaceLighting','phong','FaceColor','interp',...
        'AmbientStrength',0.7)
    material dull;
    shading interp;%faceted;
    light('Position',[0.5 0 1],'Style','infinite');
    daspect([1 1 1]);
    colorbar;
    title(text);
    colormap jet;
    if nargin<6
%         axis([0 ax 0 ay min(z) az])
        axis([min(x) max(x) min(y) max(y) 0 az])
    end
    view([-44,50])
%     view([0,0])
%     view([0 90]);
%     view([45 30]);

    if ac==1
        colorer
        colormap(terra);
%         colormap jet;
        caxis([0 az]);
    end
    if ac==2
%         cmap=rivercolor;
%         colormap(cmap);
        colormap jet;
        caxis([0 1e4]);
    end
%     caxis([0 10000]); %hacked for flacuptest12_1, 9/12/13
    if ac==7
        caxis([-1e-5 1e-5]);
        colorer;
        colormap(wave);
    end
    if ac==25
        caxis([-1e-3 1e-3]);
        colorer;
        colormap(wave);
    end
    if ac==9
        caxis([0 1]);
        colorer;
        colormap(gravelthreshold);
    end
    % Set plotting stuff
    %   if nargin>=7
    if ac==11
        caxis([-1E-4 1E-4])
    end
    if ac==13
        caxis([0 8])
    end
    if ac==15
        caxis([0 1e5])
    end
    if ac==16
        caxis([0 50])
    end
    if ac==17
        caxis([-.001 .001])
    end
%     if ac==18
%         caxis([0 5])
%     end
    %   end
    
    % Capture the frame
%     if ac==5 || ac==7
%         m(g)=getframe;
%         g=g+1;
%     else
%For coupled runs
if flacchild
    if mod(i,2) %TURN on for flacchild models
        m(g)=getframe;
        g=g+1;
    else
        continue
    end
else
    m(i)=getframe;
end
    if ac==12
        if i==2
            time=tm;
        end
        meanz(i,1)=mean(data((data(:,1)>51 & data(:,1)<65), 3));
        if i==1
            meanz(i,2)=meanz(i,1);
        else
            meanz(i,2)=mean(c2((data(:,1)>50 & data(:,1)<65), 1))*time+meanz(i-1,1);
        end
        meanz(i,3)=tm;
    end
end
if ac==12
    erate(:,1)=(meanz(2:end,2)-meanz(2:end,1))./time;
    erate(:,2)=meanz(2:end,3);
else
    erate=0;
end
% Play the movie
% movie(m);
if ac~=10
    data(:,3)=c;
else
    data=c;
end
% close the files
fclose(nfid);
fclose(tfid);
fclose(zfid);
fclose all;
