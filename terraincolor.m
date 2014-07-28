function [colorscale,shadow]=terraincolor
%generate terrain style color scale with shading to accommodate shadows in
%indexed gif

WBR = linspace(1,0.32,50) .';
WBG = linspace(1,0.19,50) .';
WBB = linspace(1,0.19,50) .';
BDgR = linspace(.32,0.11,50) .';
BDgG = linspace(.19,.31,50) .';
BDgB = linspace(.19,0.21,50) .';
DgGR = linspace(0.11,0.17,40) .';
DgGG = linspace(0.31,0.51,40) .';
DgGB = linspace(0.21,0.34,40) .';
YDbR = linspace(1,0.04,6) .';
YDbG = linspace(1,0.14,6) .';
YDbB = linspace(0.4,0.42,6) .';
GYR = linspace(0.17,1,6) .';
GYG = linspace(0.51,1,6) .';
GYB = linspace(0.34,0.4,6) .';

% WBR = linspace(1,0.32,50) .';
% WBG = linspace(1,0.19,50) .';
% WBB = linspace(1,0.19,50) .';
% BGR = linspace(.32,0,50) .';
% BGG = linspace(.19,.5,50) .';
% BGB = linspace(.19,0,50) .';
% GBlR = linspace(0,0.39,50) .'; %The green to blue adds too much blue. Try to get a "Sea Level" threshold.
% GBlG = linspace(.5,0.47,50) .';
% GBlB = linspace(0,.64,50) .';
colorscale=zeros(148,3);

% colorscale(1:50,1)=WBR;
% colorscale(1:50,2)=WBG;
% colorscale(1:50,3)=WBB;
% colorscale(50:99,1)=BGR;
% colorscale(50:99,2)=BGG;
% colorscale(50:99,3)=BGB;
% colorscale(99:148,1)=GBlR;
% colorscale(99:148,2)=GBlG;
% colorscale(99:148,3)=GBlB;

colorscale(1:50,1)=WBR;
colorscale(1:50,2)=WBG;
colorscale(1:50,3)=WBB;
colorscale(50:99,1)=BDgR;
colorscale(50:99,2)=BDgG;
colorscale(50:99,3)=BDgB;
colorscale(99:138,1)=DgGR;
colorscale(99:138,2)=DgGG;
colorscale(99:138,3)=DgGB;
colorscale(138:143,1)=GYR;
colorscale(138:143,2)=GYG;
colorscale(138:143,3)=GYB;
colorscale(143:148,1)=YDbR;
colorscale(143:148,2)=YDbG;
colorscale(143:148,3)=YDbB;


colorscale=flipud(colorscale);

%NOW THE SHADOWS
WK = linspace(1,0,20);

BKR = linspace(0.32,0,20);
BKG = linspace(0.19,0,20);
BKBl = linspace(0.19,0,20);

GKR = linspace(0,0,20);
GKG = linspace(0.5,0,20);
GKBl = linspace(0,0,20);

BlKR = linspace(.39,0,20);
BlKG = linspace(0.47,0,20);
BlKBl = linspace(0.64,0,20);

BWR = linspace(0.32,1,20);
BWG = linspace(0.19,1,20);
BWBl = linspace(0.19,1,20);

GWR = linspace(0,1,20);
GWG = linspace(0.5,1,20);
GWBl = linspace(0,1,20);

BlWR = linspace(.39,1,20);
BlWG = linspace(0.47,1,20);
BlWBl = linspace(0.64,1,20);

shadow=zeros(140,3);
shadow(1:20,1)=WK;
shadow(1:20,2)=WK;
shadow(1:20,3)=WK;

shadow(21:40,1)=BKR;
shadow(21:40,2)=BKG;
shadow(21:40,3)=BKBl;

shadow(41:60,1)=GKR;
shadow(41:60,2)=GKG;
shadow(41:60,3)=GKBl;

shadow(61:80,1)=BlKR;
shadow(61:80,2)=BlKG;
shadow(61:80,3)=BlKBl;

shadow(81:100,1)=BWR;
shadow(81:100,2)=BWG;
shadow(81:100,3)=BWBl;

shadow(101:120,1)=GWR;
shadow(101:120,2)=GWG;
shadow(101:120,3)=GWBl;

shadow(121:140,1)=BlWR;
shadow(121:140,2)=BlWG;
shadow(121:140,3)=BlWBl;



