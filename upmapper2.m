function upmapper2(basenm,basedir,step)
%This reads a flac output and interpolates the flac uplift field onto the
%child field
filenm=('child_vel979_1kmres.txt');
xyzb = creadxyzb(basenm,step);
nodes=length(xyzb);
for i=1:nodes
    xc(i) = xyzb(i,1);
    yc(i) = xyzb(i,2);
    zc(i) = xyzb(i,3);
    zci(i)=zc(i);
    bc(i) = xyzb(i,4);
end

[idf,xf,yf,zf,vf,cf] = readflac2([basedir filenm]);

vc=griddata(xf,yf,vf,xc,yc);
vcn=griddata(xf,yf,vf,xc,yc,'nearest');
vc(isnan(vc))=vcn(isnan(vc));
vc(vc<0)=0;
vc=vc';
save upmap001 vc /ascii;