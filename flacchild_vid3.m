%Flacchild_vid3: data from child side of things
function m=flacchild_vid3

fid=fopen('childdat.txt','r');
i=978;
while ~feof(fid)
    x=str2num(fgets(fid));
    y=str2num(fgets(fid));
    z=str2num(fgets(fid));
    k=str2num(fgets(fid));
    tri=delaunay(x,y);
    
    h=trisurf(tri,x,y,z,k);%,'EdgeAlpha',0.1); %edgealpha crashes matlab
    set(h,'FaceLighting','phong','FaceColor','interp',...
        'AmbientStrength',0.7)
    material dull;
    shading interp;%faceted;
    light('Position',[0.5 0 1],'Style','infinite');
    daspect([1 1 .1]);
    axis([min(x) max(x) 0 max(y) 13 3000]);
    caxis([13 2400]);
%     caxis([-4.4375 -3.3495]);
%     view([-44,50])
    view([0,90])
    colorbar;
    ac=7; colorer; colormap(wave);
% colormap jet
title([num2str(max(z)) ' m max height']);

    m(i)=getframe; 
        if i ==1
        [mrow,mcol,mwidth]=size(m(1).cdata);
    end
    if length(m(i).cdata(:,1,1))~=mrow || length(m(i).cdata(:,:,1))~=mcol
        m(i).cdata=imresize(m(i).cdata,[mrow,mcol]);
    end
i=i+1;
end