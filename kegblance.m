function [ Line_length ] = kegblance( keg_pressure,lineID,desired_tap_pressure,height )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if lineID==3/16
    resistance=3;
if nargin == 2
    Line_length=(keg_pressure-1)/resistance;
end
if nargin == 3
    Line_length=(keg_pressure-desired_tap_pressure)/resistance;
end
if nargin == 4
    Line_length=(keg_pressure-desired_tap_pressure-(height/2))/resistance;
end

end

