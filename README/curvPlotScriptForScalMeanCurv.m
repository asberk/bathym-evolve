%% Script to make log-curvature plots with zero isocontour overlays for all generalized mean curvature level surfaces

% requires *nix System
W = what('../evoByRoC/resdata');

filenames = W.mat;
savenames = cell(length(filenames), 1);

for k = 1:length(filenames)
    
    [status, savenames{k}] = system(['echo "', filenames{k}, '" | sed s/\\.mat// | sed s/\\./-/g | sed s/_/-/g | cat']);
    curvPlot(['../evoByRoC/resdata/', filenames{k}], savenames{k});
end
