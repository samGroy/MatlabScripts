limit=round(sqrt(10000.^2+10000.^2));
a=0;
rhovec=zeros(length(dat),180);
for i=1:length(dat);
    if dat(i,3)<1000 && dat(i,3)>990
            a=a+1;
        for theta=0:1:179
            rhovec(a,theta+1)=dat(i,1)*cos(deg2rad(theta))+dat(i,2)*sin(deg2rad(theta));
        end
        rhovec=round(rhovec/100)*100;
        plot(rhovec(a,:),'-k'); hold on;
    end
    a=0;
end


tri=delaunay(dat(:,1),dat(:,2));
subplot(1,2,1), trisurf(tri,dat(:,1),dat(:,2),dat(:,3));
shading interp;
subplot(1,2,2), contourf(hmatrix(:,:,2));
shading interp;
% axis([0 1000 0 180 0 max(hmatrix(:,:,2))*2])