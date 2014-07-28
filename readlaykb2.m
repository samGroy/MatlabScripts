function c=readlaykb2(basenm,ts,numg,allnodes,ac)
filesys='';
ii=3;

%read
filenm= [filesys basenm '.lay' num2str(ts-1)];
lfid=(fopen(filenm,'r'));
layfile=fscanf(lfid,'%f');
intnodes=layfile(2);
surfdat=nan(allnodes,1);

%Fetch erodibility data from the layer files.
if ac==10 % extract kb
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
elseif ac==10.1 % extract kb
        for i=1:intnodes
            if layfile(ii+6)==1
                if numg==2
                    surfdat(i)=layfile(ii+13);
                else
                    surfdat(i)=layfile(ii+12);
                end
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
elseif ac==18 % extract sediment thickness
    for i=1:intnodes
        if layfile(ii+6)==1
            surfdat(i)=layfile(ii+4);
        else
            surfdat(i)=0;
        end
        if numg==1
            ii=ii+1+layfile(ii)*7;
        elseif numg==2
            ii=ii+1+layfile(ii)*8;
        end
    end
    c=surfdat;
elseif ac==19 % bedrock/regolith flagmap
    for i=1:intnodes
        surfdat(i)=layfile(ii+6);
        if numg==1
            ii=ii+1+layfile(ii)*7;
        elseif numg==2
            ii=ii+1+layfile(ii)*8;
        end
    end
    c=surfdat;
elseif ac==20 % total exposure time
    for i=1:intnodes
        surfdat(i)=layfile(ii+3);
        if numg==1
            ii=ii+1+layfile(ii)*7;
        elseif numg==2
            ii=ii+1+layfile(ii)*8;
        end
    end
    c=surfdat;
elseif ac==21 % recent activity time
    for i=1:intnodes
        surfdat(i)=layfile(ii+2);
        if numg==1
            ii=ii+1+layfile(ii)*7;
        elseif numg==2
            ii=ii+1+layfile(ii)*8;
        end
    end
    c=log10(surfdat+1);
elseif ac==22 % creation time
    for i=1:intnodes
        surfdat(i)=layfile(ii+1);
        if numg==1
            ii=ii+1+layfile(ii)*7;
        elseif numg==2
            ii=ii+1+layfile(ii)*8;
        end
    end
    c=surfdat;
elseif ac==23 % # layers
    for i=1:intnodes
        surfdat(i)=layfile(ii);
        if numg==1
            ii=ii+1+layfile(ii)*7;
        elseif numg==2
            ii=ii+1+layfile(ii)*8;
        end
    end
    c=surfdat;
elseif ac==24 %get exposure time of bedrock. That's even if the bedrock is currently covered by sediments
    for i=1:intnodes
        if layfile(ii+6)==0
            surfdat(i)=layfile(ii+3);
            if numg==1
                ii=ii+1+layfile(ii)*7;
            elseif numg==2
                ii=ii+1+layfile(ii)*8;
            end
        elseif layfile(ii+6)==1
            if numg==1
                surfdat(i)=layfile(ii+10);
                ii=ii+1+layfile(ii)*7;
            elseif numg==2
                surfdat(i)=layfile(ii+11);
                ii=ii+1+layfile(ii)*8;
            end
        end
    end
    c=surfdat;
end

fclose(lfid);


