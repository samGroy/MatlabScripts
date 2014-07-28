function p=pmapper(basenm,pdat,opt)

%pmapper, for all your pmapping needs
%SGR--10/24/2013
%given initial CHILD run basenm and precipitation grid data pdat, pmapper
%will fit pdat to basenm and generate a precipitation map suitable for CHILD input,
%using FLOWGEN=8 and DETACHMENT_LAW=4.

%The FLOWGEN=8 option was hacked by SGR after help from NG, 10/23/2013
%DETACHMENT_LAW=4 by NG, necessary to omit all those crazy width params for erosion

%get stuff
xyzb=creadxyzb(basenm,1);
for i=1:1
    iid=fopen([basenm '.id1'],'r');
    tm=fscanf(iid,'%f',1);
    nn=fscanf(iid,'%f',1);
    ids=fscanf(iid,'%f',nn);
end
if nargin<3
    xyp=pdat; b=1;
    return
elseif strcmp(opt,'dynamic')
%     p=oro_p(basenm,2); b=0;
    p=oro_psimple(basenm,2); b=0;
%     p=2.*p;
elseif strcmp(opt,'import')
    pid=fopen(pdat,'r');
    xyp=fscanf(pid,'%f',[3,inf]);
    xyp=xyp'; b=1;
    fclose(pid);
end

%match ID codes
xyzb(:,5)=ids;

%interp stuff
if b==1
    proot=griddata(xyp(:,1),xyp(:,2),xyp(:,3),xyzb(:,1),xyzb(:,2));
    proot(isnan(proot))=0;
    %sort by descending boundary code, then by ascending index code (CHILD has odd tastes)
    [dud,index]=sortrows(xyzb,[-4 5]);
    
    %have precip data follow proper index for CHILD import
    p=proot(index,:);
end

p(p<0)=0;

save PMap.txt p -ascii
fclose(iid);