function [Xd Yd] = bilinZC(x, y, F, N, varargin)
% BILINZC - Uses bilinear interpolation to calculate approximately the zero isocontour of a given function.
%
% The function is represented by the values of the matrix F,
% corresponding to the values of the function on the points
% of the rectangular domain, which is implicitly represented by the
% vectors x and y.
%   
% F(1,1) corresponds to (x(1), y(1)); F(1,end) with (x(end), y(1));
% F(end, 1) with (x(1), y(end)); F(end, end) with (x(end), y(end))
%
% N represents the maximum number of points to interpolate in
% each domain "box", where a box is given by four mutually
% neighbouring elements of the matrix F. The line segments
% resulting from the interpolation will not, in general, be of
% uniform length; in particular, the points are evenly spaced
% in the x direction, so that the arc lengths may differ
% greatly between pairs of points. 
%
% Functionality exists for a fifth argument. This argument is a
% Boolean (represented as 1 or 0), corresponding to whether or
% not to display debugging information. Default assumption no
% debugging info. 
%
% As future update: Instead of spacing points evenly along the
% x-axis, space points according to the arc length of the zero
% contour, so that the length of each line segment defined by two
% points calculated on the isocontour is equal.
%
  
% Ensure that x and y are row vectors
  if size(x, 1)>1 & size(x,2)==1
    x = x.';
  end
  if size(y,1)>1 & size(y,2)==1
    y = y.';
  end

  % Initialize return variables
  Xd = [];
  Yd = [];
  
  % static information
  sy = length(y);
  sx = length(x);
  
  % find boxes framed by elements where a zero contour will exist
  % (i.e., we actually need the function to cross zero in order
  % for there to be a zero contour. Find these locations).
  %
  % Define a matrix whose elements are linearly independent (these
  % elements 1, 10, 100, 1000 were chosen more or less
  % arbitrarily. What's important is that each element is a
  % different magnitude than the others).
  A = [1 100; 10 1000]; % elements are L.I. by magnitude

  % example of how the convolution works for the box
  % {1,3,7,9} \in F
  % for F = [1 3 5; 7 9 13]: C_11 = 1*1 + 7*10 + 3*100 + 9*1000
  C = abs(conv2(sign(F), A, 'valid'));
  % list of element types we don't need to examine
  nonZero = [1111, 1110, 1101, 1011, 111];
  % find members of zero boxes
  LIA = ~ismember(C, nonZero);
  
  % store row and column vectors as  variables to allow 
  % references to subsets
  colvec = 1:sx-1;
  rowvec = 1:sy-1;

  % iterate over elements of F
  % ((selected columns of F containing nontrivial zero contour
  % points))
  for k = colvec(any(LIA)) 
    dx = x(k+1)-x(k);

    % set up desired domain points for current column
    xd = linspace(x(k), x(k+1), N+1);
    % handle the case when we're in the last column (to allow
    % points along far right boundary. 
    if ~(k==sx-1)
      xd = xd((2:end)-1);
    end

    % Generate Kronecker matrices
    % --------------------------------------------------------- -
    % Columns organized as Q_tl, Q_bl, Q_tr, Q_br (i.e., first 
    % column is multiplied to top-left value of function, etc.)
    % rows organized as a, c, b, d, where f(x,y)=a+bx+cy+dxy
    %
    % no longer want y(2:end) --- we want to make our choices
    % according to which boxes cross zero.  select from rowvec
    % only the elements in column k of LIA that correspond to
    % boxes where graph(F) crosses zero in non-trivial way. (the
    % plus one is used to reference the "bottom" y points)
    bases = kron([x(k+1), -x(k); -1 1],...
		 [-y(rowvec(LIA(:,k))+1).',...
		  y(rowvec(LIA(:,k))).'; 1 -1]);
    % OLD:
    % bases = kron([x(k+1), -x(k); -1 1],...
    % 		 [-y(2:end).', y((2:end)-1).'; 1 -1]);
    
    % number of elements in y matrix; = half the number of
    % rows of bases matrix; = number of boxes in current column
    nBox = sum(LIA(:, k));
    
    q = 1;
    % selected rows of F
    for j = rowvec(LIA(:,k))
      M = F(j:j+1, k:k+1);

      dy = y(j+1)-y(j);

      % will have to cleverly reference indices of Kronecker matrix:
      % q is the products of x and y {-x2y1, x2y2, x1y1, -x1y2}
      % nBox+1 is x {x2 -x2 -x1 x1}
      % nBox+1+q is y {y1 -y2 -y1 y2}
      % 2*(nBox+1) = end row is the row of ones.  {-1 1 1 -1}
      basis = bases([q, nBox+1, nBox+1+q, 2*(nBox+1)], :);
      acbd = basis*M(:);
      yd = -(acbd(1) + acbd(3)*xd)./(acbd(2) + acbd(4)*xd);
      % yd returns as +/- Inf for singular denominator
            
      % Debugging info
      if nargin>4 & varargin{1}==1
	% then debug==TRUE
	fprintf('Checking size of Xd...\n');
	fprintf(['Before: %', num2str(floor(log10(size(Xd,2)))), ...
		 'd; '], size(Xd,2));
      end

      if j == sy-1
	Xd = [Xd xd(yd >= y(j) & yd <= y(j+1))];
	Yd = [Yd yd(yd >= y(j) & yd <= y(j+1))];
      else
	Xd = [Xd xd(yd >= y(j) & yd < y(j+1))];
	Yd = [Yd yd(yd >= y(j) & yd < y(j+1))];
      end

      % Handle case where c+dx==0 --- can't divide by zero!
      % ----------------------------------------------------
      % If there are singularities, then we do the same process
      % as above, but for the inverse function:
      %               x = -(a+cx)/(b+dx)

      % (define domain in inverse space)
      yyd = linspace(y(j), y(j+1), N);
      % (check that denominator of inverse is not singular)
      if acbd(3)~=0 & any(acbd(4)*yyd~=0)
	% (obtain image points in inverse space 
	%  (i.e. domain points in regular space))
	xxd = -(acbd(1)+acbd(2)*yyd)./(acbd(3)+acbd(4)*yyd);
	
	% analogous to the last time this was done
	if k == sx-1
	  Xd = [Xd xxd(xxd >= x(k) & xxd <= x(k+1))];
	  Yd = [Yd yyd(xxd >= x(k) & xxd <= x(k+1))];
	else
	  Xd = [Xd xxd(xxd >= x(k) & xxd < x(k+1))];
	  Yd = [Yd yyd(xxd >= x(k) & xxd < x(k+1))];
	end
      end

      % increase counter by one 
      % (i.e., move focus from box q to box q+1)
      q = q+1;

      
      % Debugging info
      if nargin>4 & varargin{1}==1
	% print update on length of Xd (and thus Yd)
	fprintf(['After: %', num2str(floor(log10(size(Xd,2)))), ...
		 'd...\n'], size(Xd, 2)); 
      end
    end
  end
  
  
