function curvPlot(mfile, savename)
% assumes LStoolbox path is loaded
load(mfile);
[curvature, gradMag] = feval(@curvatureSecond, g, d_curv);
hold off;
imagesc(g.vs{2}, g.vs{1}, log(1+abs(curvature)));
colormap gray;
tmp = colormap;
colormap(tmp(end:-1:1,:));
colorbar;
hold on;
contour(g.vs{2}, g.vs{1}, d_curv, [0 0], 'b');
title('Plot of $\log~(1+|\kappa|)$ with zero isocontour overlay', 'Interpreter', 'latex', 'FontSize', 16);
xlabel('Latitude (degrees)', 'FontSize', 14);
ylabel('Longitude (degrees)', 'FontSize', 14);
set(gcf, 'PaperPosition', [0, 0, 8.5, 8.5*g.shape(1)/g.shape(2)]);
print(gcf, '-r300', '-dpng', ['~/Dropbox/academic/NSERC_Kevlahan/matlab/bathymetry/evoLS/README/figures/', savename, '.png']);
close all;