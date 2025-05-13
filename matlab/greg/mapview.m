function [] = mapview(zvar,flipcb)

% plot any map variable zvar that is either 2d, 3d, or 4d
% provided that the first two dimensions are lon and lat
% if zvar is 3d or 4d then only the first index value of dim 3 and 4 will be plotted
%
% INPUTS:
% zvar is 2d, 3d, or 4d array of the map variable with lon and lat as first two dims
% lon is vector of longitudes
% lat is vector of latitudes

if nargin==1
	flipcb = false;
end
	
% colormap('jet');
if ~flipcb
	colormap(turbo());
else
	colormap(flipud(turbo()));
end

zvar_ndims = ndims(zvar);

if zvar_ndims == 2

% check view of the data using lon -180:180 degE
figure(1)
hFig=gcf;
clf(hFig);
clear X Y Z
X = 1:size(zvar,1);
Y = 1:size(zvar,2);
Z = squeeze(zvar(:,:))';
% [C,h] = contourf(X,Y,Z,256);
% set(h,'LineColor','none')
h = pcolor(X,Y,Z);
set(h,'edgecolor','none')

elseif zvar_ndims == 3

% check view of the data using lon -180:180 degE
figure(1)
hFig=gcf;
clf(hFig);
clear X Y Z
X = 1:size(zvar,1);
Y = 1:size(zvar,2);
Z = squeeze(zvar(:,:,1))';
% [C,h] = contourf(X,Y,Z,256);
% set(h,'LineColor','none')
h = pcolor(X,Y,Z);
set(h,'edgecolor','none')

elseif zvar_ndims == 4

% check view of the data using lon -180:180 degE
figure(1)
hFig=gcf;
clf(hFig);
clear X Y Z
X = 1:size(zvar,1);
Y = 1:size(zvar,2);
Z = squeeze(zvar(:,:,1,1))';
% [C,h] = contourf(X,Y,Z,256);
% set(h,'LineColor','none')
h = pcolor(X,Y,Z);
set(h,'edgecolor','none')

end

xlabel('dim 1 index')
ylabel('dim 2 index')
h=colorbar;
set(gcf, 'PaperPosition', [0 0 8 6])   
print(gcf, [pwd '/mapview.png'], '-dpng', '-r300' ); 
