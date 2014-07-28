function result=simple_layerer(basenm,kb,ts)
% vertically continuous layers
filenm=[basenm '.lay'];
% lfid=fopen(filenm,'r');
% if lfid>0
%     fprintf('Layer file "%s" already exists.\n',filenm);
%     error('Try a different name.\n');
% end
lfid = fopen(filenm,'w');
if lfid<1
    error('Unable to create layer file.\n');
end

% Write header information
fprintf(lfid,' %.2f\n',ts);
fprintf(lfid,'%d\n',length(kb));

% Loop over nodes to write information
for j=1:length(kb)
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',1);
    
    % Basic properties: creation time, recent activity time, exposure
    % time, thickness, erodibility, and regolith/bedrock flag
    fprintf(lfid,'%.2f %.2f %.2f\n',0,0,0);
    %         if C2(layindex(i),j)==c1
    fprintf(lfid,'%.2f %f %.0f\n',10000,kb(j),0); %end
    %         if C2(layindex(i),j)==c2 || C2(layindex(i),j)==c3
    %             fprintf(lfid,'%.2f %f %.0f\n',thickness2(layindex(i),j),kb(layindex(i),j),1); end
    fprintf(lfid,'%.2f\n',10000);
end
fclose(lfid);
result=1;
