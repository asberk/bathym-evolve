function evoLSPlotSummary(matFile)
% EVOLSPLOTSUMMARY - This function generates a series of plots
% summarizing the data returned from the level set evolution,
% which is contained in matFile.
% 
% Input: 
% matFile is a path to the matlab file WITHOUT the .mat
% suffix. 
%
% Output: 
  
% This file outputs (possibly more than) one output image
% (controlled by appropriately setting the conditional statements
% below. The images display either the zero isocontour of an
% evolved level surface, obtained via the implicit
% representation of a data matrix; the level surface itself; or
% a visualization of the geometric behaviour of the level
% surface local to the zero isocontour (i.e., normal and
% tangent information).  
%
% For an example of how to run this file, see summarizeResData
% in the constCurv/ folder of the evoLS library. 
%   
  % Graphics stuff --- may be system specific
    fprintf('running administration stuff...\n');
  set(gcf, 'Visible', 'off', 'PaperUnits', 'inches');
  reslnWidth = 1440; % width of output images in pixels  
  dpi = 700;
  
  %% File administration
  savedir = 'summplots/'; % location of output files 
  % Load matFile
  load(matFile); % contains description, g, d_curv;
  [datapath rootName ext] = fileparts(matFile);
  
  %% Set up data
  fprintf('setting up data...');
  % Approximate surface normal
  [Nx, Ny, Nz] = surfnorm(d_curv);
  % rearrange domain because of how latitude and longitude were
  % input. 
  X = g.xs{2}; Y = g.xs{1}(end:-1:1, :); 
  

  %% For displaying geometric information about the level surface
  
  % We cannot display ALL geometric information about the level
  % surface, lest the visualization be horridly muddled, busy
  % and hard to read. We control the sparsity of the
  % information present (presented as a quiver vector field)
  % using the parameter "skip", refering to how many
  % rows/columns to omit from the visualization. 
  % (number of rows, columns to skip when generating sparse
  % vector field )
  skip = 8;
  % binary matrix saying if element is in sparse field
  inSparse = zeros(size(X));
  inSparse(1:skip:end, 1:skip:end) = 1;
  fprintf('success!\n');
  
  %% Generate zero isocontour. 
  simpleContour = 0;
  if simpleContour
  fprintf('generating contour...');
  [C0, h0] = contour(X, Y, d_curv, [0 0], 'b');
  print(['-r',num2str(dpi)], '-dpng', ...
	[savedir, 'isoctr0.', rootName,'.png']);
  fprintf('success!\n');
  close all;
  end
  
  
  %% This information is used in the generation of plot(s) of
  % normal vectors to the level surface that are local to the
  % zero isocontour. 
  fprintf('setting up information for quiver...');
  % define when points on d_curv are "close" to zero
  zeroThresh = 0.05; 
  % get binary matrix where 1 corresponds to that element in
  % d_curv being "close enough" to zero. 
  nearZero = abs(d_curv) < zeroThresh; 
  % Only plot the arrows that are both: far enough apart (at
  % least /skip/ elements apart), and also those that are
  % "close enough" to zero.  
  arrows = nearZero & inSparse; 
  fprintf('success!\n');
  
  % If true, the code in this conditional statement will generate
  % a 3-D quiver plot of the normal vectors to the level surface,
  % displaying only arrows that are "near" the zero isocontour.
  % (n.b.: as this plot is three-dimensional, it doesn't look
  % very good as a printed image, but may be useful in Matlab
  % if the camera can be rotated in the plot display). 
  threeDquiver = 0;
  if threeDquiver
    fprintf('generating 3D quiver...');
    % make a quiver plot of this. 
    h1 = quiver3(X(arrows), Y(arrows), d_curv(arrows), ...
		 Nx(arrows), Ny(arrows), Nz(arrows)); 
    set(h1, 'AutoScale', 'on', 'AutoScaleFactor', 0.3);
    % add the zero contour. (may want to make line width thicker
    % so the isocontour is actually visible
    hold on;
    contour(X, Y, d_curv, [0 0], 'b');
    camorbit(0, 45);
    fprintf('success!\n');
    % change paper size for a decently sized plot
    % set(gcf, 'PaperPosition',...
    % 	   [0 0 (reslnWidth/dpi) (reslnWidth/dpi)/g.shape(2)*g.shape(1  )]);
    % print plot
    fprintf('printing 3D quiver...');
    print(['-r',num2str(dpi)], '-dpng',...
	  [savedir, 'ZeroIsoCtrNormalArrows.', rootName, '.png']);
    hold off; close all;

    fprintf('success!\n');
  end

  % if true, the code in this conditional statement generates a
  % 2-D quiver plot of the vector field that is normal to the
  % level surface, only displaying arrows that are "near" the
  % zero isocontour of the level surface. 
  twoDquiver = 1;
  if twoDquiver
    fprintf('generating 2D quiver...');
    h1 = quiver(X(arrows), Y(arrows), Nx(arrows), Ny(arrows));
    set(h1, 'AutoScale', 'on', 'AutoScaleFactor', 0.8);
    % add the zero contour. (may want to make line width thicker
    % so the isocontour is actually visible
    % hold on;
    % contour(X, Y, d_curv, [0 0], 'b');
    % camorbit(0, 45);
    fprintf('success!\n');
    % change paper size for a decently sized plot
    % set(gcf, 'PaperPosition',...
    % 	   [0 0 (reslnWidth/dpi) (reslnWidth/dpi)/g.shape(2)*g.shape(1  )]);
    % print plot
    fprintf('printing 2D quiver...');
    print('-r900', '-dpng',...
	  [savedir, 'ZeroIsoProjdNormls.', rootName, '.png']);
    hold off; close all;
    fprintf('success!\n');
  end

  fprintf('\nexiting function call.\n\n');  
  
