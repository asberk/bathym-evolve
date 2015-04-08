% This script uses the default Matlab function contour to generate a zero
% contour of the bathymetry data contained in images/etopo1Large.tif.
% It uses latitude & longitude data contained, by default, in the .tif file

% run this code from the README directory

geoA = imread('../images/etopo1Large.tif');
load '../images/geoRLarge.mat';
X = geoR.Lonlim(1)+geoR.DeltaLon/2:geoR.DeltaLon:geoR.Lonlim(2)-geoR.DeltaLon/2;
Y = geoR.Latlim(2)-geoR.DeltaLat/2:geoR.DeltaLat:geoR.Latlim(1)-geoR.DeltaLat;

contour(X,Y, geoA, [ 0 0 ], 'b');

pbaspect([size(geoA,2), size(geoA,1), 1]);