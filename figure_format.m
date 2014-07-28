figure
X=0; Y=X;
set(gcf,'PaperUnits','centimeters')
xSize = 8.9; ySize = xSize+1;
xLeft = (21.59-xSize)/2; yTop = (27.94-ySize)/2;
set(gcf,'PaperPosition',[xLeft yTop xSize ySize]);
set(gcf,'Position',[X Y xSize*50 ySize*50]);
x = .15; y = .25;
width = 1-(x+.05); height = 1-(y+.05);
axes1=axes('position',[x y width height]);
color=colorbar('SouthOutside');
set(color,'position',[x 0.09 width 0.04]);

% Geology 2-column width = 8.9 cm
% Geology page width = 18.5 cm