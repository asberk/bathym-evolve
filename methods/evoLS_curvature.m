%% data = evoLS_curvature(data,g,tMax, bVal, accuracy,evoByRoC)
%  
% Update summary: 2013-09-03
% update mostly includes better commenting, and a little
% cleaning. 
% Update summary: 2013-08-23
% update included removal of RoCFlow in favour of incorporation
% of termRadiusOfCurvature. evoByRoC and useRoCDependent
% parameters were retained. 
% Update summary: 2013-08-07
% update included the initial writing of RoCFlow, as well as
% its incorporation in the body of evoLS_curvature. Argument
% splitFlow (originally from curvatureStarDemo) was modified to
% evoByRoC, and useTimeDependent was changed to
% useRoCDependent. 
% Previous updates included adding additional arguments and
% re-structuring of the base algorithm. 

function data = evoLS_curvature(data, g, tMax, bVal, accuracy, evoByRoC)
% evoLS_curvature: method used to evolve a given level surface by
% mean curvature (by default) or "generalized" mean curvature
% (i.e., subject to some scaling multiplier). 
%
%   [ data ] = evoLS_curvature(data, g, tMax, bVal, accuracy, evoByRoC)
%
% Parameters:
%   accuracy: Controls the order of approximations.  Note that the spatial
%   approximation is always second order.
%
%                  'low'         Use odeCFL1.
%                  'medium'      Use odeCFL2 (default).
%                  'high'        Use odeCFL3.
%
%   evoByRoC: Boolean.  Use the radius of curvature-dependent
%   version of the flow field, which adjusts the relative b value
%   coefficient in proportion to the local radius of curvature (in
%   fact, by the local curvature, itself).
%   g: struct.  The domain on which the surface is to be
%   evolved. All grids should be passed _after_ passing them
%   through the function processGrid from the LStoolbox. 
%   data: double.  A matrix object that represents the height
%   of the initial level surface at each point in space (thus
%   assuming a regular [Cartesian] grid cell structure). 
%   bVal: scalar. this value is used as the constant multiplier
%   to the curvature term in the level set equation. Good
%   results have been obtained with this value in the range of
%   0.02--0.04. 
%   tMax: scalar. This represents the duration of the level set
%   evolution. Good results have been obtained for tMax in the
%   range of 1--2.5. 
%
% Output Parameters:
%   data: Implicit surface function at tMax.
%

% Original file copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% curvatureStarDemo and LStoolbox written 
% by Ian Mitchell, 2/13/04;
% 
% This file almost entirely modified from the original 
% by Aaron Berk (aberk@math.toronto.edu), 2013. 
  
% %---------------------------------------------------------
% Note: need to have LStoolbox "Kernel" folder in MATLAB's
% path, but inside this function is not the right place to
% add it. 
% The following should be in any script file that uses this
% function: 
% % Make sure we can see the kernel m-files.
% run('../addPathToKernel');

%-----------------------------------------------------------
% Setting parameter values
% Duration parameter.
if(nargin < 3) 
  tMax = 1.75;
end

% Curvature speed parameter.
if(nargin < 4)
  bValue = 0.02;
else
  bValue = bVal;
end

if(nargin < 5)
  accuracy = 'medium';
end

% Use the time dependent motion?
if(nargin < 6)
  useRoCDependent = 0;
else
  useRoCDependent = evoByRoC;
end

%-----------------------------------------------------------
% Integration parameters.

% Can change this vector if desired. 
plot_points = 0:tMax/10:tMax; % length = 16
%tMax = max(plot_points);                  % End time.
t0 = 0;                                   % Start time.

% How close (relative) do we need to get to tMax to be considered
% finished?
small = 100 * eps;

%-----------------------------------------------------------
% Set up functions used for motion by mean curvature

schemeData.grid = g;
schemeData.curvatureFunc = @curvatureSecond;

if(useRoCDependent)
  % Time dependent flow field
  % Use modified termCurvature function to increase speed (and
  % hopefully accuracy, somehow?)
  schemeFunc = @termRadiusOfCurvature;
  schemeData.b = bValue;
  schemeData.tSwitch = 0.5 * tMax;
else
  schemeFunc = @termCurvature;
  % Time independent flow field is constant.
  schemeData.b = bValue;
end

%-----------------------------------------------------------
% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.9, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  integratorFunc = @odeCFL1;
 case 'medium'
  integratorFunc = @odeCFL2;
 case 'high'
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

%-----------------------------------------------------------
% Initialize time stepping (variables are poorly named)
plot_count = 1;

if(plot_points(1) == t0)
  plot_count = plot_count + 1;
end

%-----------------------------------------------------------
% Loop until tMax (up to allowed roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, plot_points(plot_count) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

  plot_count = plot_count + 1;
  
end

endTime = cputime;
fprintf('Total execution time %g seconds\n', endTime - startTime);
end
