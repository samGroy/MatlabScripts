ind=(rand(size(x))>0.9);
xs=x(ind); ys=y(ind); zs=z(ind);

[xu,xind]=unique(xs);
[yu,yind]=unique(ys);
zmat=nan(numel(yind),numel(xind));
zind=sub2ind([yind xind],size(zmat));

xmat=repmat(xu',numel(yind),1);
ymat=repmat(yu,1,numel(xind));
zmat(ind)=zs;

figure
subplot(121)
scatter(xs,ys,10,zs,'filled')
axis image

subplot(122)
pcolor(xmat,ymat,zmat)
shading flat
axis image