function [stmatrix,stcov]=stminterpolate(stmin,t,pntCrd,varargin)
%STMINTERPOLATE   Spatio-temporal interpolation of a space time matrix.
%  STMATRIX=STMINTERPOLATE(STMIN,T,PNTCRD) interpolates the space-time
%  matrix structure STMIN, and returns a new space-time matrix (matrix,
%  not structure) for the decimal years in T and points in PNTCRD. PNTCRD
%  must contain the latitude, longitude (and height) of the points in the
%  output space-time matrix STMATRIX. STMIN is a space-time matrix
%  structure, but STMATRIX is a plain [numPoints x numEpochs x ndim]
%  matrix.
%
%  [STMATRIX,STCOV]=STMINTERPOLATE(...) also returns the covariance
%  matrix STCOV.
%
%  [...]=STMINTERPOLATE(...,'option',value,...) allows to specify options
%  for the interpolation
%
%    'spatial'    spatial interpolation method ['triangulation']
%    'temporal'   temporal interpolation method (any method supported by
%                 interp1), default 'linear'
%    'pntCrdType' Point coordinate type ['deg/m','km/m'] (default 'deg/m')
%
%  See also interp1, delaunayTriangulation and scatteredInterpolant.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:   8 Oct 2020 by Hans van der Marel
% Modified: .. Oct 2020 by Hans van der Marel
%            - changes something

% Check input arguments and process options

if nargin < 3
   error('This function expects at least three input arguments.')
end

opt.spatial='triangulate';
opt.temporal='linear';
opt.pntCrdType='deg/m';
for k=1:2:numel(varargin)
    opt.(varargin{k})=varargin{k+1};
end
    
% IMPORTANT Temporal interpolation is not yet supported, so the output times must
% match the input matrix.
%
% We can use Matlabs interp1 with any of the methods 
%
% Question: what to do, first spatial, or first temporal interpolation....

numEpochs=size(t);
if any( abs(t-stmin.epochDyear) > 1e-6 )
   error('This function does not yet support temporal interpolation.')
end


% Convert latitude and longitude to km for interpolation

switch lower(opt.pntCrdType)
   case 'deg/m'
      % convert latitude/longitude into local topocentric coordinates
      [pntNeu,plh0] = plh2neusp(pntCrd);          % deg/deg/m -> m/m/m  OUT  
      [tmpNeu] = plh2neusp(stmin.pntCrd,plh0);  % deg/deg/m -> m/m/m  IN
      pntNeu(:,1:2)=pntNeu(:,1:2)./1000;          % m/m/m -> km/km/m
      tmpNeu(:,1:2)=tmpNeu(:,1:2)./1000;          % m/m/m -> km/km/m
   case 'km/m'
      % coordinates are already in the right units, this is exceptional,
      % and only happens for simulations, just copy
      pntNeu=pntCrd;
      tmpNeu=stmin.pntCrd;
   otherwise
      error('unknown pntCrdType option')        
end
      
% Interpolate to the points at hand (only in space, we assume the
% epochs are the same (to be fixed!)

switch opt.spatial
    case 'triangulate'

        % The following code fails when point are not inside the convex hull
        % or on it's border (ti and bc contain Inf entries). We replace it by
        % Matlabs scatteredInterpolant, but this does not allow error
        % propagation.

        % % Create a delaunay Triangulation of the scattered points 
        %
        % DT = delaunayTriangulation(tmpNeu(:,1:2));
        % % Find the triangle that encloses each query point using the pointLocation 
        % % method. In the code below, ti contains the IDs of the enclosing triangles 
        % % and bc contains the barycentric coordinates associated with each triangle.
        % [ti,bc] = pointLocation(DT,pntNeu(:,1:2));
        % % Find the point id's of the vertices surrounding the triangle
        % triIdx=DT(ti,:);
        % % Calculate the sum of the weighted values of V(x,y) using the dot product.
        % for k=1:numEpochs
        %  for l=1:3
        %     v=stmin.obsData(:,k,l);
        %     stmatrix(:,k,l)=dot(bc,v(triIdx),2);
        %  end
        % end

        % The previous code fails when points are outside the triangles. This is
        % unfortunate, because we need the points numbers in order to do the 
        % error propagations.
        % 
        % Possible solutions
        % 1. when ti and bc have a NaN (outside triangle), use only the nearest
        %    neighbor instead of interpolating between three neighbors
        % 2. use Matlabs scatteredInterpolant function, which does linear
        %    extrapolation by default. This is the best way to go for the actual
        %    interpolation, but unfortunately, this class does not return the 
        %    triangulation, and therefore, we cannot do error propagation.
        % Way to go. Use 2. for interpolation, use 1. for error propagation. For
        % error propagation it matters if we do interpolation or extrapolation.

        % Build the interpolant object (using the height as temporary value)
        F=scatteredInterpolant(tmpNeu(:,1:2),tmpNeu(:,3));
        % Do the actual interpolation 
        for k=1:numEpochs
           for l=1:3
              F.Values=apriori.obsData(:,k,l);
              stmatrix(:,k,l)=F(pntNeu(:,1:2));
           end
        end
        
    otherwise
        error('unknown interpolation method')
end
        
end