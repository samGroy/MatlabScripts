%attempt to interpolate surfaces from baselay, faster rendering compared to
%point data
giggity=baselay';
coh=nan(length(giggity),max_num_layers);
xx=xc'; yy=yc';
snarf=zeros(length(giggity),1);
laysurf=nan(length(giggity),max_num_layers);
for i=1:length(giggity)
    snarf(i)=sum(~isnan(giggity(i,:)));
    if snarf(i)==1
        laysurf(i,max_num_layers)=giggity(i,1);
        coh(i,3)=layer_cohesion(1,i);
    end
    if snarf(i)==2
        laysurf(i,2)=giggity(i,1);
        laysurf(i,3)=giggity(i,2);
        coh(i,2)=layer_cohesion(1,i);
        coh(i,3)=layer_cohesion(2,i);
    end
    if snarf(i)==3
        laysurf(i,1)=giggity(i,1);
        laysurf(i,2)=giggity(i,2);
        laysurf(i,3)=giggity(i,3);
        coh(i,1)=layer_cohesion(1,i);
        coh(i,2)=layer_cohesion(2,i);
        coh(i,3)=layer_cohesion(3,i);
    end
end

tri=delaunay(xx,yy);
trisurf(tri,xx,yy,laysurf(:,1),coh(:,1));
hold on;
trisurf(tri,xx,yy,laysurf(:,2),coh(:,2));
hold on;
trisurf(tri,xx,yy,laysurf(:,3),coh(:,3));
hold off;
