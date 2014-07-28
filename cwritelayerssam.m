% Create layer file and open it for writing (first make sure it doesn't
% already exist)

basenm='test01';
ts=0;

filenm= [ basenm '.lay' num2str(ts) ]; %replace test01 with basenm once I get my sh!t together
lfid=fopen(filenm,'r');
% if lfid>0
%     fprintf('Layer file "%s" already exists.\n',filenm);
%     error('Try a different name.\n');
% end
lfid = fopen(filenm,'w');
if lfid<1
    error('Unable to create layer file.\n');
end

ALLUVIUM_FLAG=1;

% Write header information
fprintf(lfid,' %.2f\n',ts);
fprintf(lfid,'%d\n',length(bc));

% Loop over nodes to write information
for j=1:length(xc)
    
    % Write number of layers at this node
    fprintf(lfid,' %.0f\n',nlayers(j));

    % For each layer at this node, write information
    for i=1:nlayers(j)
        
        % Basic properties: creation time, recent activity time, exposure
        % time, thickness, erodibility, and regolith/bedrock flag
        fprintf(lfid,'%.2f %.2f %.2f\n',0,0,0);
        fprintf(lfid,'%.2f %f %.0f\n',layer_thickness(i,j),kb(i,j),0);
        
%         % Grain size information
%         if numg==1
%             fprintf(lfid,'%.2f\n',laydat(j,i,1));   % this line gives 'old' format: if one size, just write thickness
%         elseif numg>1
%             fprintf(lfid,'%.2f ',laydat(j,i,1)-sum(laydat(j,i,5:(3+numg))));  % this line gives 'old' format: if more than one, write thickness of size 1 as total thick minus sum of sizes 2-N
%             for k=2:numg-1
%                 fprintf(lfid,'%f ',laydat(j,i,4+(k-1)) );
%             end
%             fprintf(lfid,'%f\n',laydat(j,i,4+(numg-1)) );
%         end
        
    end
    
end