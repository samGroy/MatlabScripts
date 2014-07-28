limit=round(sqrt(10000.^2+10000.^2));
rhovec=[0:99]'; rhovec=rhovec*(limit/100);
% rhoint=100;
% thetaint=10;
hmatrix=zeros(length(rhovec),180); hmatrix2=hmatrix;
for i=1:length(dat);
    for theta=0:1:179
        rho=abs(dat(i,1)*cos(deg2rad(theta))+dat(i,2)*sin(deg2rad(theta)));
        if rho<=limit && rho >=0
            n=abs(rho-rhovec);
            [e,d]=min(n);
            hmatrix(d,theta+1)=hmatrix(d,theta+1)+dat(i,3);
            hmatrix2(d,theta+1)=hmatrix2(d,theta+1)+1;
        end
    end
end
c=hmatrix./hmatrix2;
tri=delaunay(dat(:,1),dat(:,2));
subplot(1,2,1), trisurf(tri,dat(:,1),dat(:,2),dat(:,3));
shading interp;
view([0 90]);
x=[0:179];
subplot(1,2,2), surf(x,rhovec,c);
view([0 90]);
% shading interp;
% axis([0 1000 0 180 0 max(hmatrix(:,:,2))*2])