function c=readlaykb3(basenm,ts,numg,allnodes,ac)
filesys='';
ii=3;

%read
filenm= [filesys basenm '.lay' num2str(ts-1)];
lfid=(fopen(filenm,'r'));
layfile=fscanf(lfid,'%f');
intnodes=layfile(2);
surfdat=nan(allnodes,1);

%Fetch erodibility data from the layer files.
for i=1:intnodes
    if layfile(ii+6)==1
        surfdat(i)=1;
    else
        surfdat(i)=layfile(ii+5);
    end
    if numg==1
        ii=ii+1+layfile(ii)*7;
    elseif numg==2
        ii=ii+1+layfile(ii)*8;
    end
end
c=log10(surfdat);

fclose(lfid);


