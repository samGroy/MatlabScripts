function [knick_d,knick_h]=knickpt_tracker(distance,height, height_orig)
%knickpt_tracker: the latest attempt at finding knickpoints
% Uses cstrmproseries_vid to get upstream distance and height data.
% Options:
% 1) measures slope of the data over a wide search area, and the greatest
% slope equals the knickpoint position.
% 2) compares first time step profile to the current time step, finds
% position where certain difference threshold exists (if knickpoint is 10m,
% call it 1m)

% This option is for slope
% slope=zeros(1,length(height));
% for i=5:1:length(height)-5
%     slope(i)=(height(i-4)-height(i+5))/(distance(i-4)-distance(i+5));
% end
% 
% [val,ind]=min(slope);
% 
% knick_d=distance(ind);
% knick_h=height(ind);
    
drop=5;

% This option is for height difference thresholding
dh=height-height_orig';
for i=1:length(dh)
    if dh(i)<= -drop 
        ind=i; 
        break 
    end
    if i==length(dh)
        ind=length(dh);
    end
end


if dh(1)<=-drop
    knick_d=0;
    knick_h=height_orig(ind)-drop;
elseif ind==length(dh)
    knick_d=distance(ind);
    knick_h=height_orig(ind-1)-drop;
else
    slope=(dh(ind-1)-dh(ind))/(distance(ind-1)-distance(ind));
    proper_d=(-drop-dh(ind-1))/slope;
    proper_d=proper_d+distance(ind-1);
    knick_d=proper_d;
    knick_h=height_orig(ind)-drop;
end
