function p=pmapperFLAC(yc,zc)

%pmapper, for all your pmapping needs
%SGR--10/24/2013
%given initial CHILD run basenm and precipitation grid data pdat, pmapper
%will fit pdat to basenm and generate a precipitation map suitable for CHILD input,
%using FLOWGEN=8 and DETACHMENT_LAW=4.

%The FLOWGEN=8 option was hacked by SGR after help from NG, 10/23/2013
%DETACHMENT_LAW=4 by NG, necessary to omit all those crazy width params for erosion

ridgept=yc(zc==max(zc)); % Try defining a ridgepoint for now...think of something later...10/2013
ridgept=mean(ridgept);
p(yc>=ridgept)=.05;
p(yc<ridgept)=.001;

save PMap.txt p -ascii
% fclose(iid);