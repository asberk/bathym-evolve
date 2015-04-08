%% This is an example script of how to use bilinZC to obtain
% the zero contour using bilinear interpolation of a given
% surface function, which is represented as a regular matrix. 

% load file from constCurv/resdata as example. 
% contains objects: d_curv, g, description 
load('../../parentPath.mat');
load([parentPath, ...
      'evoLS/constCurv/resdata/mapL_curvatureLS.tMax2.3.b0.04.mat']);
% import bilinZC
addpath(genpath([parentPath, 'evoLS/methods']));

% set domain. (Lat and Lon were input incorrectly, so domain
% is weird by default and we have to 'fix' it.)
x = g.vs{2};
y = g.vs{1}(end:-1:1, :);

% set maximum number of points to create in each grid cell box.
% note that resulting grid points will not, in general, be
% evenly spaced. We could change this, but it's more work that
% may not be necessary. 
N = 5; 

% execute bilinear interpolation of zero contour.
[Xd, Yd] = bilinZC(x, y, d_curv, N);

% Note that MarkerSize can be adjusted to 5 to match the
% linewidth of the contour (5 =/= 5px)
figure(1);
% MATLAB's interpolation
[~, hc] = contour(x, y, d_curv, [0 0], 'r');
% fix aspect ratio
pbaspect([g.shape(2), g.shape(1) 1]);
hold on;
% plot bilinear interpolation
plot(Xd, Yd, '.b', 'MarkerSize', 7);
% change layering of points and lines
uistack(hc, 'top')
hold off;
% UX
legend('Bilinear interpolation approximation', 'Matlab''s contour function', 'Location', 'North');
title('Bilinear interpolation approximation to zero isocontour of level surface', 'FontSize', 16);
xlabel('Latitude (degrees)');
ylabel('Longitude (degrees)');

makeFigure = 0;
if makeFigure
    set(gcf, 'PaperUnits', 'inches', 'PaperPosition', [0 0 8.5 8.5*g.shape(1)/g.shape(2)]);
    print(gcf, '-r400', '-dpng', './../README/figures/bilin-interp-approx.png');
    close all;
end

