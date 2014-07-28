for i=1:length(coords(:,1))
    dx=coords(i,1)-coords(:,1);
    dy=coords(i,2)-coords(:,2);
    dist=sqrt(dx.^2+dy.^2);
    ind=find(dist<=5e3);
    sigmaT(i)=std(tilt(ind));
    sigmaAR(i)=std(aspect_ratio(ind));
end