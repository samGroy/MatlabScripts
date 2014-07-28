limit=round(sqrt(10000.^2+10000.^2));
sublimit=-limit;
rhovec=[sublimit:limit]';
hmatrix=zeros(length(rhovec),180);
for i=1:10:length(dat);
    for theta=0:1:179
        rho=dat(i,1)*cos(deg2rad(theta))+dat(i,2)*sin(deg2rad(theta));
        rho=round(rho);
        if rho<=limit && rho >=sublimit
            d=find(rhovec==rho);
            hmatrix(d,theta+1)=hmatrix(d,theta+1)+dat(i,3);
        end
    end
end

tri=delaunay(dat(:,1),dat(:,2));
subplot(1,2,1), trisurf(tri,dat(:,1),dat(:,2),dat(:,3));
shading interp;
subplot(1,2,2), image(hmatrix);
% shading interp;
% axis([0 1000 0 180 0 max(hmatrix(:,:,2))*2])