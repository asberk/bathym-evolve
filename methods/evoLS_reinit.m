function data = AtlSignedDistFn(data, g, accuracy)
% reinitDemo: demonstrate the reinitialization equation.
%
% [ data, g, data0 ] = reinitDemo(initialType, accuracy, displayType)
%
% Demonstrates the application of the reinitialization equation to turn one
%   of several different dynamic surface functions into signed distance
%   functions.  While it works in any dimension, it is hard to visualize
%   the difference in dimensions higher than two.
%  
% This function was originally designed as a script file, so most of the
%   options can only be modified in the file.
%
% For example, edit the file to change the grid dimension, boundary conditions,
%   flow field parameters, etc.
%
% Parameters:
%
%   initialType  String to specify which initial dynamic implicit surface.
%                  'circle'      An off center circle/sphere (default).
%                  'star'        A star-shaped interface (not 1D).
%   accuracy     Controls the order of approximations.
%                  'low'         Use odeCFL1 and upwindFirstFirst (default).
%                  'medium'      Use odeCFL2 and upwindFirstENO2.
%                  'high'        Use odeCFL3 and upwindFirstENO3.
%                  'veryHigh'    Use odeCFL3 and upwindFirstWENO5.
%   displayType  String to specify how to display results.
%                  The specific string depends on the grid dimension;
%                  look at the helper visualizeLevelSet to see the options
%                  (optional, default depends on grid dimension).
%
%   data         Implicit surface function at t_max.
%   g            Grid structure on which data was computed.
%   data0        Implicit surface function at t_0.

% Copyright 2004 Ian M. Mitchell (mitchell@cs.ubc.ca).
% This software is used, copied and distributed under the licensing 
%   agreement contained in the file LICENSE in the top directory of 
%   the distribution.
%
% Ian Mitchell, 2/13/04

%---------------------------------------------------------------------------
% You will see many executable lines that are commented out.
%   These are included to show some of the options available; modify
%   the commenting to modify the behavior.
  
%---------------------------------------------------------------------------
% Make sure we can see the kernel m-files.
%run('../addPathToKernel');

%---------------------------------------------------------------------------
% Integration parameters.
tMax = 1.0;                  % End time.
tPlot = 0.1;                 % Period at which plot should be produced.
t0 = 0;                      % Start time.
%singleStep = 0;              % Plot at each timestep (overrides tPlot).

% How close (relative) do we need to get to tMax to be considered finished?
small = 100 * eps;

%---------------------------------------------------------------------------
if(nargin < 3)
  accuracy = 'medium';
end

% Set up time approximation scheme.
integratorOptions = odeCFLset('factorCFL', 0.5, 'stats', 'on');

% Choose approximations at appropriate level of accuracy.
switch(accuracy)
 case 'low'
  derivFunc = @upwindFirstFirst;
  integratorFunc = @odeCFL1;
 case 'medium'
  derivFunc = @upwindFirstENO2;
  integratorFunc = @odeCFL2;
 case 'high'
  derivFunc = @upwindFirstENO3;
  integratorFunc = @odeCFL3;
 case 'veryHigh'
  derivFunc = @upwindFirstWENO5;
  integratorFunc = @odeCFL3;
 otherwise
  error('Unknown accuracy level %s', accuracy);
end

%if(singleStep)
%   integratorOptions = odeCFLset(integratorOptions, 'singleStep', 'on');
% end

% Set up spatial approximation scheme.
schemeFunc = @termReinit;
schemeData.grid = g;
schemeData.derivFunc = derivFunc;
schemeData.initial = data;
% Use the subcell fix by default.
schemeData.subcell_fix_order = 1;

%---------------------------------------------------------------------------
% Loop until tMax (subject to a little roundoff).
tNow = t0;
startTime = cputime;
while(tMax - tNow > small * tMax)

  % Reshape data array into column vector for ode solver call.
  y0 = data(:);

  % How far to step?
  tSpan = [ tNow, min(tMax, tNow + tPlot) ];
  
  % Take a timestep.
  [ t y ] = feval(integratorFunc, schemeFunc, tSpan, y0,...
                  integratorOptions, schemeData);
  tNow = t(end);

  % Get back the correctly shaped data array
  data = reshape(y, g.shape);

end

endTime = cputime;
fprintf('\nTotal execution time %g seconds\n', endTime - startTime);
