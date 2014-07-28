% Local correlogram map scale vid
% Given a bunch of points, each with numerous aspect ratio and azimuth
% values at discreet scales, this script incrementally takes the data,
% plots, advances the scale, and repeats. 
% Also assumes you have the base map that you ned up and ready to go.
v=aspect_ratio; %What you're plotting
freezeColors;
hold on; colormap jet;
caxis([0 1]);
tri=delaunay(coords(:,1),coords(:,2));
for i=1:length(tilt(1,:))
    trisurf(tri,coords(:,1),coords(:,2),coords(:,3)+100,v(:,i));
    shading interp; material dull;
%     scatter3(coords(:,1),coords(:,2),coords(:,3)+50,100,tilt(:,i),'filled')
    m(i)=getframe;
end

% You can then save the video as avi.
