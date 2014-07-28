function kegbalance( keg_pressure,lineID,desired_tap_pressure,height )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if lineID==3/16
    resistance=3;
end
if desired_tap_pressure>1
    Line_length=desired_tap_pressure;
    desired_tap_pressure=(-Line_length*resistance)+keg_pressure-(height/2);
    fprintf('\nKeg Pressure: %f PSI\nTap Pressure: %f PSI\nLine Length: %f feet\nLine Resistance: %f PSI/foot\nHeight to Tap: %f feet\nCalibrated to 40 F\n\n', keg_pressure,desired_tap_pressure,Line_length,resistance,height);
elseif nargin == 2
    Line_length=(keg_pressure-1)/resistance;
    fprintf('\nKeg Pressure: %f PSI\nTap Pressure: 1 PSI\nLine Length: %f feet\nLine Resistance: %f PSI/foot\nCalibrated to 40 F\n\n', keg_pressure,Line_length,resistance);
elseif nargin == 3
    Line_length=(keg_pressure-desired_tap_pressure)/resistance;
    fprintf('\nKeg Pressure: %f PSI\nTap Pressure: %f PSI\nLine Length: %f feet\nLine Resistance: %f PSI/foot\nCalibrated to 40 F\n\n', keg_pressure,desired_tap_pressure,Line_length,resistance);
elseif nargin == 4
    Line_length=(keg_pressure-desired_tap_pressure-(height/2))/resistance;
    fprintf('\nKeg Pressure: %f PSI\nTap Pressure: %f PSI\nLine Length: %f feet\nLine Resistance: %f PSI/foot\nHeight to Tap: %f feet\nCalibrated to 40 F\n\n', keg_pressure,desired_tap_pressure,Line_length,resistance,height);
end
end

