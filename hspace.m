function hspace(dat,dimension)
dat(:,1)=dat(:,1)-min(dat(:,1)); dat(:,2)=dat(:,2)-min(dat(:,2));
limit(1)=max(dat(:,1))/2;
limit(2)=max(dat(:,2))/2;
limit=min(limit);
l=length(dat(:,1));
if dimension==1, l=180; end
if dimension==2
    if l<=12000
        l=100;
    else
        l=200; end
end

rhovec=[0:l-1]'; rhovec=rhovec*(2*limit/l);%change back to 0:199 and /200 for model domains%rhovec=[0:199]'; rhovec=rhovec*(2*limit/200);
hmatrix=zeros(length(rhovec),180); hmatrix2=hmatrix;
for i=1:length(dat);
    for theta=0:1:179
%         fprintf('angle %f \n',theta);
        rho=limit+((dat(i,1)-limit)*cos(deg2rad(theta))+(dat(i,2)-limit)*sin(deg2rad(theta)));
        if rho<=2*limit && rho >=0
            n=abs(rho-rhovec);
            [e,d]=min(n);
            hmatrix(d,theta+1)=hmatrix(d,theta+1)+dat(i,3);
            hmatrix2(d,theta+1)=hmatrix2(d,theta+1)+1;
        end
    end
end
if dimension==2
    tri=delaunay(dat(:,1),dat(:,2));
    subplot(1,2,1), trisurf(tri,dat(:,1),dat(:,2),dat(:,3));
    shading interp;
    axis([min(dat(:,1)) max(dat(:,1)) min(dat(:,2)) max(dat(:,2)) min(dat(:,3)) max(dat(:,3))])
    view([0 90]);
else
    subplot(1,2,1), plot(dat(:,1),dat(:,2),'.'); %comment out for anything but line data
    axis([0 max(dat(:,1)) 0 max(dat(:,2))]);
end
x=[0:179];
if dimension==1
    c=hmatrix2;
    subplot(1,2,2), surf(x,rhovec,c); %switch hmatrix2 with c for anything but line data
else
    c=hmatrix./hmatrix2;
    subplot(1,2,2), surf(x,rhovec,c);
end
% colormap(flipud(gray));
colormap(gray);
shading flat;
view([0 90]);
axis([0 max(x) 0 limit*2]);