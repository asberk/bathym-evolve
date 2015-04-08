%% N.B. This file is DEPRECATED. 
% Instead, use radiusOfCurvatureLSevo.m with "regular" evolution
% by curvature ( either by setting tSwitch to a value that is
% greater than tMax or by setting the scaling function \psi\beta
% (scalingFn) to the unity function.

%% mapLRunEvoLS will evolve the level surface of the large map
% found in geotiffFN (below) using a reinitialization of
% bathymetry data and then using mean curavture evolution to
% transform the surface.  

% This file is the same as NCSALargeCurvatureLSEvo.m except for
% the fact that this does not have the ability to load data (which
% is superfluous because the file exits after the run, anyway.

%-----------------------------

%% This function will consolidate other methods, and use 'medium' accuracy
%  curvature level set evolution to acquire a smoother version of the
%  northern coast of South America and the Central America coastal region.

%---------------------------------------------------------------------------
% Make sure we can access functions from LS package
load('../../parentPath.mat');
addpath(genpath([parentPath, 'LStoolbox/Kernel']));

%--------------------------------------------------
% add custom methods to path
addpath(genpath([parentPath, 'evoLS/methods']));

%----------------------------------------------------
% Script parameters

% size of map to evolve:
mapSize = 'Small'; % could be Large. 

% relative location of save directory in which to save output
savedir = 'resdata/';
% This is the savefile prefix for the series of images generated
% by evolving the surface of the (_mapSize_) (map) of Northern
% South America extending up to Floridian Coast using curvature
% (L)evel (S)et methods that evolve the surface at a rate
% proportional to the local curvature (inversely proportional to
% the radius of curvature).  To this prefix will be appended the
% duration of surface evolution, the curvature parameter value,
% and appropriate file extension (.mat, .png)
% 
% if doSave then this is the name of the file to save; if
% ~doSave, then this is the base name of the file to load. 
savefileprefix = ['map', num2str(mapSize(1)),'_curvatureLS']; 


% ------------------------------------------------------------
% Import and scale GeoTIFF

% relative location of image file
geotiffFN = ['../images/etopo1', mapSize, '.tif'];

try 
  fprintf('Trying to read in geotiff file using geotiffread...');
  [geoA, geoR] = geotiffread(geotiffFN);
  fprintf('success!\n');
catch ME
  fprintf('failed!\nIt''s likely that this version of MATLAB does not have geotiffread (see explanation below...)\n\n');
  display(ME);

  fprintf(['\n\nTrying alternative read in; requires geoR',...
	   mapSize, '.mat...\n']);
  if exist(['../images/geoR', mapSize, '.mat'])==2
    fprintf(['    found geoR', mapSize, '.mat\n']);

    fprintf(['    reading geoR', mapSize, '.mat\n'])
    load(['geoR', mapSize, '.mat']);
    
    fprintf('    reading TIFF image using imread...');
    geoA = imread(geotiffFN);
    fprintf('Success!\n\n');
  else % geoR file not found. 
    error(['geoR', mapSize, '.mat not found. Ensure that geoR', ...
	   mapSize, '.mat is located in ../images and try again']);
  end
end

% Convert geoA to proper format
geoA = double(geoA);
d_init = geoA./max(abs(geoA(:))); % in range (-1, 1);
		     % zero contour remains unshifted
				  
% ------------------------------------------------------------
% Make grid

g.dim = 2;
g.min = [geoR.Latlim(1) - geoR.DeltaLat/2;
	 geoR.Lonlim(1) + geoR.DeltaLon/2];
g.max = [geoR.Latlim(2) + geoR.DeltaLat/2;
	 geoR.Lonlim(2) - geoR.DeltaLon/2];
g.N = size(geoA).';
g.bdry = @addGhostExtrapolate;

g = processGrid(g);

clear geoA; % geoA not need anymore 

% -----------------------------------------------------------
%% Evolve Level Surfaces
fprintf('beginning reinitialization...\n');
d_reinit = evoLS_reinit(d_init, g); % medium accuracy:
				    % 516.74 s on 949 x
				    % 1562 = 1.48 M pixels
clear d_init; % d_init not needed anymore
fprintf('reinitialization complete...\n');

% set parameter values for curvature evolution
bVal = [0.01, 0.02, 0.04, 1];
tMax = [1, 1.5, 1.9, 2.3];
parMat = combvec(tMax, bVal);
clear bVal tMax; % not needed anymore

% check multiparm (if so, expect long run?) 
nRuns = size(parMat, 2);

for k = 1:nRuns
  fprintf('beginning curvature evolution for tMax=%5.3f, bVal=%4.3f...\n', parMat(1, k), parMat(2,k));  
  d_curv = evoLS_curvature(d_reinit, g, parMat(1,k), parMat(2,k)); % medium accuracy:
     % 354.45 s on 949 x
     % 1562 = 1.48 M
     % pixels
  fprintf('curvature evolution %02d complete...\n', k);
  
  % save result
  description = ['This file was created using the script mapLRunEvoLS.m using medium accuracy for both reinitialization and level set evolution by curvature. The duration of this run was ', num2str(parMat(1,k)), ', with a bValue of ', num2str(parMat(2,k)),'. (Note that d_init can be obtained by rescaling the original image found in the associated tiff file: ', geotiffFN, ' .'];
  % define name of savefile. 
  savefile = [savedir,savefileprefix, '.tMax',num2str(parMat(1,k)), '.b', num2str(parMat(2,k)), '.mat'];
  save(savefile, 'description', 'g', 'd_curv');
  
  fprintf('Successfully saved results for run %02d...\n\n', k);
  clear d_curv; %remove before next run. 
end

fprintf('Runs finished...\n\n');

% Stuff for remote runs to quit the matlab once run is completed. 
fprintf('\nexiting script...\n');
W = what;
if ~isempty(strfind(W.path, '/1/home/berkas/'));
  exit;
end
