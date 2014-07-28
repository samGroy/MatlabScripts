function m=flacchild_vid(basenm,output,opt_no_vid)
%flacchild vid: used to create videos from the saved flacdat files from a
%coupled simulation. basenm is the name of the GIF video you want to make,
%output is 1) elevation, 2) uplift rate, 3) cohesion. If added, opt_no_vid
%doesn't save the GIF of the movie.

%Set parameters
keepgoing=1; counter=1;
figure('Color',[1 1 1]);
zmin=0; zmax=1000;%zmin=0; zmax=1500; %These values may differ between runs
if output==1, cmin=0; cmax=900; end
if output==2, cmin=-.05; cmax=.05; end
if output==3, cmin=6; cmax=7.4471; end
fullcolorscale=zeros(288,3);

while keepgoing == 1
    filesys='';
    filenm= ['child_vel' num2str(counter) '.txt' ];
    nfid=fopen(filenm,'r');
    if nfid<=0, keepgoing=0; break,end
    nn = fscanf(nfid,'%d',1);
    idxyzwc=rot90(fscanf(nfid,'%f',[6,nn]));
    tri=delaunay(idxyzwc(:,2),idxyzwc(:,3));
    if output==1
        h=trisurf(tri,idxyzwc(:,2),idxyzwc(:,3),idxyzwc(:,4)); %elev is color
        [cmap,shadow]=terraincolor; %get colormap
    elseif output==2
        h=trisurf(tri,idxyzwc(:,2),idxyzwc(:,3),idxyzwc(:,4),idxyzwc(:,5)); %uplift/w is color
        cmap=jet;
    elseif output==3
        h=trisurf(tri,idxyzwc(:,2),idxyzwc(:,3),idxyzwc(:,4),log10(idxyzwc(:,6))); %log10(cohesion) is color
        cmap=jet;
    end
    % update counter, set frame stuff, get frame
    view([30 40]);
    daspect([1 1 .02]);
    %     axis([min(idxyzwc(:,2)) max(idxyzwc(:,2)) min(idxyzwc(:,3)) max(idxyzwc(:,3)) zmin zmax]);
    axis([-1000 max(idxyzwc(:,2)) min(idxyzwc(:,3)) max(idxyzwc(:,3)) zmin zmax]);
    caxis([cmin cmax]);
    colormap(cmap);
    set(h,'FaceLighting','phong','FaceColor','interp',...
        'AmbientStrength',0.7)
    material dull;
    shading interp;
    light('Position',[1 1 .25],'Style','infinite');
    %fullcolorscale(1:148,1:3)=cmap;
    %fullcolorscale(149:288,1:3)=shadow;%add the white background and black mesh/axes to the colormap!
    m(counter)=getframe;
    %Uncomment if framing size trouble occurs....
    if counter ==1
        [mrow,mcol,mwidth]=size(m(1).cdata);
    end
    if length(m(counter).cdata(:,1,1))~=mrow || length(m(counter).cdata(:,:,1))~=mcol
        m(counter).cdata=imresize(m(counter).cdata,[mrow,mcol]);
    end
    %GIF function shut down for now, can't handle shadows...
%     if nargin==2
%         % Name GIF, make the GIF vid
%         if output==1, gifname=([filesys basenm '_elev.gif']); end
%         if output==2, gifname=([filesys basenm '_uplift.gif']); end
%         if output==3, gifname=([filesys basenm '_cohesion.gif']); end
%         
%         im=frame2im(m(counter));
%         imind = rgb2ind(im,cmap);
%         if counter==1
%             imwrite(imind,cmap,gifname,'gif','DelayTime',0.07,'Loopcount',inf);
%         else
%             imwrite(imind,cmap,gifname,'gif','DelayTime',0.07,'WriteMode','append');
%         end
%     end
    fclose(nfid);
    counter=counter+1;
end
if nargin == 3
    if opt_no_vid == 2
        if output==1, obj=VideoWriter([basenm '_elev']); end
        if output==2, obj=VideoWriter([basenm '_uplift']); end
        if output==3, obj=VideoWriter([basenm '_cohesion']); end
        obj.Quality=100;
        obj.FrameRate=20;
        open(obj);
        writeVideo(obj,m);
        close(obj);
    end
end

    
    
