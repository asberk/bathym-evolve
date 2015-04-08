%% radiusOfCurvatureLSevo.m

%% This function will act as a framework to consolidate other
% methods in the implementation of a level set method.  By
% default, the level set method is written to use evolution by
% radius of curvature. By default, this script runs method using
% 'medium' accuracy, and the "large" etopo1 geotiff file. 
%
% Note that many parts of this file can be modified to add or
% remove functionality. 

% Last modified: 3 September 2013. 

%-----------------------------------------------------------
% Add paths. 
%
% This script was written with the use of rsync in mind so that
% the entire evoLS folder can be mirrored across different
% systems. Consequently, it is expected that there is a .mat file
% containing a single string object which is located in the same
% directory as the evoLS folder. The string object contains the
% absolute path to the "Kernel" folder containing the methods of
% the LStoolbox. Note that this string object can also be used as
% an absolute reference to the location of the evoLS folder.
load('../../parentPath.mat');
addpath(genpath([parentPath, 'LStoolbox/Kernel']));

%--------------------------------------------
% add path to custom methods
addpath(genpath([parentPath, 'evoLS/methods']));

%-------------------------------------------
% Script parameters

% size of map to evolve:
mapSize = 'Large'; % {'Small', 'Large'}. 

% relative location of save directory in which to save output
savedir = 'resdata/';

% The savefile prefix is used for the series of images generated
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
savefileprefix = ['map', mapSize(1),'_RoCLS']; 


% ----------------------------------------------------------
% Import and scale GeoTIFF

geotiffFN = [parentPath, 'evoLS/images/etopo1', mapSize, '.tif'];
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
    load(['../images/geoR', mapSize, '.mat']);
    
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
    
%-----------------------------------------------------------
% Make grid
    
g.dim = 2;
g.min = [geoR.Lonlim(1) + geoR.DeltaLon/2;
	 geoR.Latlim(1) - geoR.DeltaLat/2];
g.max = [geoR.Lonlim(2) - geoR.DeltaLon/2;
	 geoR.Latlim(2) + geoR.DeltaLat/2];
g.N = size(geoA).';
g.bdry = @addGhostExtrapolate;

g = processGrid(g);

clear geoA; % geoA not needed anymore
    
%-----------------------------------------------------------
% Evolve Level Surfaces proportional to local curvature

fprintf('Attempting reinitialization...\n');

% medium accuracy: 516.74 s on 949 x 1562 = 1.48 M pixels
d_reinit = evoLS_reinit(d_init, g); 
clear d_init; 
fprintf('Reinitialization complete.\n');

fprintf('Starting level set evolution(s) by RoC...\n');
fprintf('--parametrizing runs...\n');
% note that bVal and tMax can be vectors. 
bVal = 0.04;
tMax = 2.3;
parMat = combvec(tMax, bVal);
clear bVal tMax; % not needed anymore. 

% check multiparm (if so, expect long run?) 
nRuns = size(parMat, 2);

fprintf('--total number of runs is %02d\n', nRuns);

for k = 1:nRuns
  fprintf(['--beginning curvature evolution %02d of %02d;\n' ...
	   '----tMax=%5.3f, bVal=%4.3f...\n'], ...
	  k, nRuns, parMat(1, k), parMat(2,k));  
  % medium accuracy: 354.45 s on 949 x 1562 = 1.48 M pixels
  d_curv = evoLS_curvature(d_reinit, g, parMat(1,k), ...
			   parMat(2,k), 'low', 1);
  fprintf('---Complete!\n');
  
  fprintf('--Saving result...');
  
  description = ['An attempt at evolving etopo1', mapSize, ...
		 ' using evoLS_curvature --- a function' ...
		 ' that attempts to evolve the surface' ...
		 ' according to local radius of curvature' ...
		 ' (if the surface is ''too curved'',' ...
		 ' then the surface evolves more quickly;' ...
		 ' if surface curvature is within' ...
		 ' appropriate bounds, then evolution is' ...
		 ' locally suppressed). Evolution was run' ...
		 ' via radiusOfCurvatureLSevo.m. The simulation' ...
		 ' duration of the run was', num2str(parMat(1,k)),...
		 ' with a bMltplr of', num2str(parMat(2,k)), ...
		 '. Note that d_init can be obtained by' ...
		 ' rescaling the original image found in the' ...
		 ' associated tiff file: ', geotiffFN, ').']; 
  
  % define name of savefile. 
  savefile = [savedir, savefileprefix, ...
	      '.tMax', num2str(parMat(1,k)), ...
	      '.b', num2str(parMat(2,k)), '.mat'];
  save(savefile, 'description', 'g', 'd_curv');
  
  clear description d_curv; % remove before next run. 
  
  fprintf('success!\n');
  
end

fprintf('Level set evolutions complete.\n\n');

% Stuff for remote runs to quit the matlab once run is completed. 
fprintf('\nexiting script...\n');
W = what;
if ~isempty(strfind(W.path, '/1/home/berkas/'));
  exit;
end

