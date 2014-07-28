filesys='';
numg=2;
duration=100;
basenm='texture02_1';
allnodes=6384;
resdat=nan(allnodes,duration);
seddat=nan(allnodes,1);
oseddat=seddat;
age=seddat;
oage=seddat;
for ts=1:duration
    ii=3;
    filenm= [filesys basenm '.lay' num2str(ts-1)];
    lfid=(fopen(filenm,'r'));
    layfile=fscanf(lfid,'%f');
    intnodes=layfile(2);
    for i=1:intnodes
        seddat(i)=layfile(ii+6);
        age(i)=layfile(ii+3);
        if ts>1
            if seddat(i)==0 && oseddat(i)==1
                resdat(i,ts)=oage(i);
            end
        end
        oseddat(i)=seddat(i);
        oage(i)=age(i);
        if numg==1
            ii=ii+1+layfile(ii)*7;
        elseif numg==2
            ii=ii+1+layfile(ii)*8;
        end
    end
end
residence_intervals=nanmean(resdat,2);

