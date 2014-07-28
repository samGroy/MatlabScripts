function faultsynch(childrunbasenm,xf,yf,zf,pf,cf)
% faultsynch: synch damage zones in CHILD to the weak zones that develop in
% strain softening FLAC

% Where does it go: replaces old opt_dynamic_erody == 2 at line 281, varres

% Steps:
% 1) read in plastic strain data and cohesion data
% 2) find plastic strain peaks that occur in strain weakened zones (if
% coh<3e7) using the findpeaks tool, seems to work ok
% 3) Check to see if the peaks alrady exist, if so keep the old builds
% 4) If new zones have formed, add a damage zone to the list
%   a) x coordinate is always x axis middle of model
%   b) y coordinate is coordinate of the peak
%   c) z coordinate is coordinate at peak
%   d) strike is always due east for convergent models
%   e) Dip is + wehn north of y midpoint, - when south, always 30
%   f) thickness is function of peak intensity maybe?
%   g) Optional build in later: gauge degree of weakening in CHILD models to
% the degree of cohesion drop in FLAC (can't do this right now)
% 5) build damage zones in synthlay using the current surface.

% Necessary inputs:
% flac data: x,y,z,c,ps
% if max(ps)>0
% search ps data along x centerline, get y coordinates, cohesion values, ps
% values
% findpeaks


if max(pf)>0
    if strainswitch==0
        fprintf('Permanent crustal strain initiated!');
        strainswitch=1;
        faultfile=zeros(1,7);
        c1=1e5; c2=2e5; c3=1e7; c4=3e7;
        thick=yc_res*2;
        strike=0; dip=30;
    end
    yline=yf(xf<100 && cf<c4); pline=pf(xf<100 && cf<c4); zline=zf(xf<100 && cf<c4);
    [pks,locs]=findpeaks(pline);
    fprintf('%d faults located\n',length(pks));
    ycoords=yline(locs); zcoords=zline(locs);
    for i=1:length(pks)
        if ycoords(i)>Ly/2
            faultfile(i+1,1)=-dip;
        else
            faultfile(i+1,1)=dip;
        end
        faultfile(i+1,2)=strike;
        faultfile(i+1,3)=thick;
        faultfile(i+1,4)=0;
        faultfile(i+1,5)=Lx/2;
        faultfile(i+1,6)=ycoords(i);
        faultfile(i+1,7)=zcoords(i);
    end
    synthlay_grad7( childrunbasenm, faultfile, Dt, 1, 1, 1, 1, 1, 1) % Should automatically spit out a .lay file for use
end
