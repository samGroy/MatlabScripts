function [c,flagmatrix]=covermap(basenm,ts,numg,allnodes,flagmatrix)
filesys='';
ii=3;

% do creadxyzb stuff, i.e. get stuff from .nodes and .z files.
% read in .nodes, .z.

filenm= [filesys basenm '.lay' num2str(ts-1)];
lfid=(fopen(filenm,'r'));
layfile=fscanf(lfid,'%f');
intnodes=layfile(2);
surfdat=nan(allnodes,1);

%     Fetch erodibility data from the layer files.
for i=1:intnodes
    surfdat(i)=layfile(ii+6);
    if numg==1
        ii=ii+1+layfile(ii)*7;
    elseif numg==2
        ii=ii+1+layfile(ii)*8;
    end
end
flagmatrix(:,ts)=surfdat;
c=mean(flagmatrix,2).*100;
fclose(lfid);


