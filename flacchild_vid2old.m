function m=flacchild_vid2(basenm,output,opt_no_vid)
%flacchild vid: used to create videos from the saved flacdat files from a
%coupled simulation. basenm is the name of the GIF video you want to make,
%output is 1) elevation, 2) uplift rate, 3) cohesion. If added, opt_no_vid
%doesn't save the GIF of the movie.

%Set parameters
keepgoing=1; counter=1;
figure('Color',[1 1 1]);
zmin=0; zmax=1000; %These values may differ between runs
if output==1, cmin=0; cmax=450; end
if output==2, cmin=-.1; cmax=.1; end
if output==3, cmin=6; cmax=7.4471; end

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
    elseif output==2
        h=trisurf(tri,idxyzwc(:,2),idxyzwc(:,3),idxyzwc(:,4),idxyzwc(:,5)); %uplift/w is color
    elseif output==3
        h=trisurf(tri,idxyzwc(:,2),idxyzwc(:,3),idxyzwc(:,4),log10(idxyzwc(:,6))); %log10(cohesion) is color
    end
    % update counter, set frame stuff, get frame
    view([30 70]);
    axis([min(idxyzwc(:,2)) max(idxyzwc(:,2)) min(idxyzwc(:,3)) max(idxyzwc(:,3)) zmin zmax]);
    caxis([cmin cmax]);
    cmap=colormap; %get colormap
    cmap(65,:)=[1 1 1]; cmap(66,:)=[0 0 0]; %add the white background and black mesh/axes to the colormap!
    m(counter)=getframe;
    if nargin==2
        % Name GIF, make the GIF vid
        if output==1, gifname=([filesys basenm '_elev.gif']); end
        if output==2, gifname=([filesys basenm '_uplift.gif']); end
        if output==3, gifname=([filesys basenm '_cohesion.gif']); end
        
        im=frame2im(m(counter));
        imind = rgb2ind(im,cmap);
        if counter==1
            imwrite(imind,cmap,gifname,'gif','DelayTime',0.07,'Loopcount',inf);
        else
            imwrite(imind,cmap,gifname,'gif','DelayTime',0.07,'WriteMode','append');
        end
    end
    counter=counter+1;
    fclose(nfid);
end

