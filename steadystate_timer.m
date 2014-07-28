% steadystate_timer: finds the limits of the apparent fault exposure at time 0,
% gets the average elevation of all points in that area, and tracks the
% elevation change as the model reaches "steady-state".

% Needs: basenm to get the faultfile specs for dip and width data. Assumes
% the fault strikes due N.

% Get the faultfile stuff
basenm='sskp06';
iteration=4;
total_ts=21;
cd(['C:\child_n\ChildExercises\paramsweep01\' basenm]);
load faultfile.mat

dip=abs(faultfile(2,1,iteration));
true_width=faultfile(2,3,iteration);
x_pos=faultfile(2,5,iteration);

% calculate the x range you'll be using to collect elevation data
x_west=x_pos;
apparent_width=1/sin(deg2rad(dip))*(5*true_width);
x_east=apparent_width+x_west;

% Get elevation data
zavg=zeros(total_ts,1);
for i=1:total_ts
    xyz=creadxyz([basenm '_' num2str(iteration)],i);
    z=xyz(xyz(:,1)>x_west & xyz(:,1)<x_east,3);
    zavg(i)=mean(z);
end