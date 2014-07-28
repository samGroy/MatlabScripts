hmatrix=zeros(1000,180,2);
for i=1:length(dat);
    for theta=0:1:180
        rho=dat(i,1)*cos(deg2rad(theta))+dat(i,2)*sin(deg2rad(theta));
        rho=round(rho);
        if rho<=1000
            hmatrix(rho,theta,1)=hough(rho,theta,1)+1;
            hmatrix(rho,theta,2)=hmatrix(rho,theta,2)+dat(i,3);
        end
    end
end
hmatrix(:,:,2)=hmatrix(:,:,2)./hmatrix(:,:,1);