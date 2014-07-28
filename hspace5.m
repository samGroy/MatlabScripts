% hough space is a good method for determining directional dependence and position for
% large scale features, no good for splitting signal by length scale.
% Preference for large scale features is because of cumulate frequency,
% greater abundance of points representing a single large feature will show stronger
% than scattered small features.
% Yields frequency of angle and radius from the center of the domain.
% Parameter can be anything spit out by CHILD.
limit=11147;
rhovec=[0:99]'; rhovec=rhovec*(2*limit/100);%change back to 0:199 and /200 for model domains%rhovec=[0:199]'; rhovec=rhovec*(2*limit/200);
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
c=hmatrix./hmatrix2;
tri=delaunay(dat(:,1),dat(:,2));
% subplot(1,2,1), trisurf(tri,dat(:,1),dat(:,2),dat(:,3)); %comment out
% when using line data
subplot(1,2,1), plot(dat(:,1),dat(:,2),'.'); %comment out for anything but line data
shading interp;
view([0 90]);
x=[0:179];
subplot(1,2,2), surf(x,rhovec,hmatrix2); %switch hmatrix2 with c for anything but line data
colormap(flipud(gray));
shading flat;
view([0 90]);
% shading interp;
% axis([0 1000 0 180 0 max(hmatrix(:,:,2))*2])