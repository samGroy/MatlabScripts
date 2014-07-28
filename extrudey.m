%load synthesized zx cohesion planes and extrude them into y
load dip50_3.dat
dip50=dip50_3;
kb=dip50(:,3);
kb(kb>=5E6)=1E7;
kb(kb<5E6)=1E4;
dip50(:,4)=kb;
dip50(:,3)=dip50(:,2);
dip50(:,2)=zeros;

slice=length(dip50);
nyslice=4;
y=slice/nyslice;
for i=1:slice
    if dip50(i,3)==0
        dip50(i,5)=1;
    else
        dip50(i,5)=0;
    end
end
%dip(slice+1:nyslice*slice,:)=zeros;
% repmat(dip50,nyslice*slice,5);
% for i=slice:slice:nyslice*slice
%     dip50(i:i+slice,2)=y;
% %     dip(i,1)=dip(1,1);
% %     dip(i:i+slice,3)=dip(1:slice,3);
% %     dip(i:i+slice,4)=dip(1:slice,4);
% %     dip(i:i+slice,5)=dip(1:slice,5);
%     y=i+1;
% end

    