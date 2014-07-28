function [v,c]=readlaykb(basenm,n,numg,az)
% basenm='fractureset10';
filesys='';
% n=51;
i=1;
ii=3;
% numg=1;
% totalnodes=9604;%81;%2401;%9801;;  %case-by-case, find way to automatically determine node count
%x=zeros(totalnodes); y=zeros(totalnodes); z=zeros(totalnodes);




for ts=1:n
    % do creadxyzb stuff, i.e. get stuff from .nodes and .z files.
    % read in .nodes, .z.
    filenm= [filesys basenm '.nodes' ];
    nfid=fopen(filenm,'r');
    filenm= [filesys basenm '.z' ];
    zfid=fopen(filenm,'r');
    for iii=1:ts
        tm = fscanf(nfid,'%f',1);
        if iii==ts,fprintf('Reading data for time %f...',tm);end
        nn = fscanf(nfid,'%d',1);
        j=fscanf(nfid,'%f',[4,nn]);
        tm = fscanf(zfid,'%f',1);
        nn = fscanf(zfid,'%d',1);
        zed=fscanf(zfid,'%f',[1,nn]);
    end
    
%     m = [ rot90(j(1,:),3) rot90(j(2,:),3) rot90(zed,3) rot90(j(4,:),3) ];
%     fclose(nfid);
%     fclose(zfid);
%     for i=1:length(m)
%         if m(i,4)==0
%             x(i,1)=m(i,1); y(i,1)=m(i,2); z(i,1)=m(i,3);
%         end
%     end
    
    
    filenm= [filesys basenm '.lay' num2str(ts-1)];
    lfid=(fopen(filenm,'r'));
    layfile=fscanf(lfid,'%f');
    totalnodes=layfile(2);
    kbsurfdat=nan(totalnodes,n);
    x=rot90(j(1,1:totalnodes),3);
    y=rot90(j(2,1:totalnodes),3);
    z=rot90(zed(1:totalnodes),3);
    
%     Fetch erodibility data from the layer files.
    if numg==1
        for i=1:totalnodes
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+1+layfile(ii)*7;
        end
    elseif numg==2
        for i=1:totalnodes
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+1+layfile(ii)*8;
        end
    end
    
    fprintf('rocks are fun!\n');
    fclose(lfid);
    i=1;
    ii=3;
    %  plot
    kb=log10(kbsurfdat(:,ts));
    tri =delaunay(x,y);
    h=trisurf(tri,x,y,z,kb);
    daspect([1 1 .1]);
%     axis([0 ax 0 ay -1000 az]);
    axis([min(x) max(x) min(y) max(y) 0 az]);
        set(h,'FaceLighting','phong','FaceColor','interp',...
        'AmbientStrength',0.7)
    material dull;
    shading interp;
    light('Position',[0.5 0 1],'Style','infinite');
%     view([0 90]);
    view([-45 60]);
    daspect([1 1 2]);
%     view([-45,60]);
    colorbar;
    colormap jet;
%     caxis([-4.4 -2.5]);
%     shading flat;
    v(ts)=getframe;
end

% movie(v);
c(:,1)=x;
c(:,2)=y;
c(:,3)=kbsurfdat(:,ts);
%dlmwrite('dip50.kb',kbsurfdat);
fprintf('now its done!\n');