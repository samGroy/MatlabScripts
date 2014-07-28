function v=readlaykb(basenm,n,totalnodes,numg,ax,ay,az)
% basenm='fractureset10';
filesys='';
% n=51;
i=1;
ii=3;
% numg=1;
% totalnodes=9604;%81;%2401;%9801;;  %case-by-case, find way to automatically determine node count
kbsurfdat=nan(totalnodes,n);
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
    
    x=rot90(j(1,1:totalnodes),3);
    y=rot90(j(2,1:totalnodes),3);
    z=rot90(zed(1:totalnodes),3);
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
    nodecount=layfile(2);
    for i=1:totalnodes
        if layfile(ii)==1
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+8;
            if numg==2, ii=ii+1; end
            continue
        end
        if layfile(ii)==2
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+15;
            if numg==2, ii=ii+2; end
            continue
        end
        if layfile(ii)==3
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+22;
            if numg==2, ii=ii+3; end
            continue
        end
        if layfile(ii)==4
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+29;
            if numg==2, ii=ii+4; end
            continue
        end
        if layfile(ii)==5
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+36;
            if numg==2, ii=ii+5; end
            continue
        end
        if layfile(ii)==6
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+43;
            if numg==2, ii=ii+6; end
            continue
        end
        if layfile(ii)==7
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+50;
            if numg==2, ii=ii+7; end
            continue
        end
        if layfile(ii)==8
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+57;
            if numg==2, ii=ii+8; end
            continue
        end
        if layfile(ii)==9
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+64;
            if numg==2, ii=ii+9; end
            continue
        end
        if layfile(ii)==10
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+71;
            if numg==2, ii=ii+10; end
            continue
        end
        if layfile(ii)==11
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+78;
            if numg==2, ii=ii+11; end
            continue
        end
        if layfile(ii)==12
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+85;
            if numg==2, ii=ii+12; end
            continue
        end
        if layfile(ii)==13
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+92;
            if numg==2, ii=ii+13; end
            continue
        end
        if layfile(ii)==14
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+99;
            if numg==2, ii=ii+14; end
            continue
        end
        if layfile(ii)==15
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+106;
            if numg==2, ii=ii+15; end
            continue
        end
        if layfile(ii)==16
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+113;
            if numg==2, ii=ii+16; end
            continue
        end
        if layfile(ii)==17
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+120;
            if numg==2, ii=ii+17; end
            continue
        end
        if layfile(ii)==18
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+127;
            if numg==2, ii=ii+18; end
            continue
        end
        if layfile(ii)==19
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+134;
            if numg==2, ii=ii+19; end
            continue
        end
        if layfile(ii)==20
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+141;
            if numg==2, ii=ii+20; end
            continue
        end
        if layfile(ii)==21
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+148;
            if numg==2, ii=ii+21; end
            continue
        end
        if layfile(ii)==22
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+155;
            if numg==2, ii=ii+22; end
            continue
        end
        if layfile(ii)==23
            kbsurfdat(i,ts)=layfile(ii+5);
            ii=ii+162;
            if numg==2, ii=ii+23; end
            continue
        end
    end
    fprintf('rocks are fun!\n');
    fclose(lfid);
    i=1;
    ii=3;
    %  plot
    kb=log10(kbsurfdat(:,ts));
    tri=delaunay(x,y);
    trisurf(tri,x,y,z,kb);
    axis([0 ax 0 ay 0 az]);
    view([-45,60]);
    colorbar;
    colormap jet;
    v(ts)=getframe;
end

movie(v);

%dlmwrite('dip50.kb',kbsurfdat);
fprintf('now its done!\n');