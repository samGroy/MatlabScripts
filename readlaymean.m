function [c,flagmatrix]=readlaymean(basenm,ts,ac,numg,allnodes,flagmatrix)
filesys='';
ii=3;


filenm= [filesys basenm '.lay' num2str(ts-1)];
lfid=(fopen(filenm,'r'));
layfile=fscanf(lfid,'%f');
intnodes=layfile(2);
surfdat=nan(allnodes,1);

%     Fetch erodibility data from the layer files.
for i=1:intnodes
    if ac==19
        surfdat(i)=layfile(ii+6)*100;
    elseif ac==18
        if layfile(ii+6)==1
            surfdat(i)=layfile(ii+4);
        else
            surfdat(i)=0;
        end
    elseif ac==21 % recent activity time
        surfdat(i)=layfile(ii+2);
    end
    if numg==1
        ii=ii+1+layfile(ii)*7;
    elseif numg==2
        ii=ii+1+layfile(ii)*8;
    end
end
flagmatrix(:,ts)=surfdat;
c=mean(flagmatrix,2);
fclose(lfid);


