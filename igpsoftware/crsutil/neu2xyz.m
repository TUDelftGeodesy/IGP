function [xyz,R] = neu2xyz(neu,plh,mode)
%NEU2XYZ Local coordinates (North,East,Up) to cartesian coordinates (dX,dY,dZ).
%   DXYZ=NEU2XYZ(NEU,PLH) converts the M-by-3 matrix NEU with local coordinates 
%   (North, East, Up) into a M-by-3 matrix DXYZ with cartesian coordinate 
%   differences (dX, dY and dZ) in a ECEF system. PLH contains the geographic 
%   coordinates, i.e. latitude and longitude (height is not used) in radians, 
%   of the center of the local system. PLH is either a M-by-2 or M-by-3 matrix, 
%   or a 1-by-2 or 1-by-3 vector. In case PLH is a vector the same center point
%   is used for all coordinate triplets.
%
%   This function has an optional third argument MODE. Default for MODE is 'p'
%   ('plh'). Other MODE's are 'n' ('normal'), 'r' and 'xr':
%
%   DXYZ=NEU2XYZ(NEU,N,'n') uses the normal vector N instead of geographic 
%   coordinates PLH of the center point. N is a M-by-3 matrix or 1-by-3 vector 
%   vector pointing in the up direction (N does not have to be a unit vector).
%
%   DXYZ=NEU2XYZ(NEU,XR,'r') uses the cartesian coordinates XR of the center 
%   point instead of geographic coordinates PLH. XR is a M-by-3 matrix or 
%   1-by-3 vector. 
%
%   XYZ=NEU2XYZ(NEU,XR,'xr') outputs cartesian coordinates XYZ (instead of 
%   coordinate differences DXYZ), and uses the center point XR as with the 
%   previous option. XYZ is computed inside the function from XR.
%
%   [DXYZ,R]=NEU2XYZ(...) also outputs the rotation matrix R, with
%   DXYZ'=R*NEU'. R is a 3-by-3 matrix in case PLH, N or XR is a vector. In 
%   case PLH, N, or XR is a matrix, R is 3-by-3-by-M. The rotation matrix
%   can be used for instance to convert the covariance matrix (when available
%   and needed).
%
%   Examples:
%       dxyz=neu2xyz(neu,[52*pi/180 4*pi/180])
%       dxyz=neu2xyz(neu,plh)
%       dxyz=neu2xyz(neu,n,'n')
%       dxyz=neu2xyz(neu,xr,'r')
%       xyz=neu2xyz(neu,xr,'xr')
%       [xyz,R]=xyz2neu(neu,xr,'xr')
%
%   Please note that the following two examples do not give the same results
%   because the up-direction is different (geodetic versus astronomic)
%
%       xr=plh2xyz(plh) 
%       dxyz=neu2xyz(neu,xr,'r')    (same result as dxyz=neu2xyz(neu,plh))
%       dxyz=neu2xyz(neu,xr,'n')
%
%   See also XYZ2NEU, XYZ2PLH, XYZ2ZAS, ZAS2XYZ and PLH2XYZ.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995,2013.

%   Created:     7 May 1995 by Hans van der Marel
%   Modified:   14 Jun 2013 by Hans van der Marel
%                - faster implementation reducing use of trigonometry 
%                - option to use normal vector N instead of phi/lambda
%                - added the options to use XR
%                - updated description
%               12 March 2014 by Hans van der Marel
%                - added the option to output the rotation matrix
%                - minor corrections to description

% Input argument checking

if nargin < 2 | nargin > 3
  error('Must be called by two or three arguments.');
end
if nargin == 2
  mode='plh';
end

[m,n]=size(neu);
if n~=3 & m==3
  neu=neu';
end
[m1,n1]=size(neu);
[m2,n2]=size(plh);
if n2~=3 & n2~=2 & ( m2==3 | m2==2 )
  plh=plh'; 
  [m2,n2]=size(plh);
end
if m1~=m2 & m2~=1
   error(['PLH, N or XR must be same length as NEU or a vector']);
end

% Sort out the various modes

switch lower(mode)
   case {'normal','n'}
      % plh contains normal vector -> normalize to unit vector
      n=plh ./ sqrt( plh(:,1).*plh(:,1) + plh(:,2).*plh(:,2) + plh(:,3).*plh(:,3) );
   case {'plh','p'}
      % plh contains latutude and longitude -> compute unit normal vector
      n= [ cos(plh(:,1)).*cos(plh(:,2)) cos(plh(:,1)).*sin(plh(:,2))  sin(plh(:,1)) ];
   case {'r','xr'}
      % plh contains XR -> compute latitude and longitude, then unit normal vector
      xr=plh;
      plh=xyz2plh(plh);
      n= [ cos(plh(:,1)).*cos(plh(:,2)) cos(plh(:,1)).*sin(plh(:,2))  sin(plh(:,1)) ];
   case {'legacy'}
      % Old way of doing things (using less efficient and flexible trigonometry)
      xyz = [ -sin(plh(:,1)).*cos(plh(:,2)).*neu(:,1) - sin(plh(:,2)).*neu(:,2) + cos(plh(:,1)).*cos(plh(:,2)).*neu(:,3) ...
              -sin(plh(:,1)).*sin(plh(:,2)).*neu(:,1) + cos(plh(:,2)).*neu(:,2) + cos(plh(:,1)).*sin(plh(:,2)).*neu(:,3) ...
               cos(plh(:,1))               .*neu(:,1)                           + sin(plh(:,1))               .*neu(:,3) ];
      if n~=3 & m==3, xyz=xyz';, end
      return
   otherwise
      error('unknown mode')
end

% Do the actual computation 

% xyz = [  -n(:,1).*n(:,3)./cphi.*neu(:,1)  - n(:,2)./cphi.*neu(:,2) + n(:,1).*neu(:,3)   ...
%          -n(:,2).* n(:,3)./cphi.*neu(:,1) + n(:,1)./cphi.*neu(:,2) + n(:,2).*neu(:,3)   ...
%           cphi.*neu(:,1)                                           + n(:,3).*neu(:,3)    ]

cphi= sqrt(1-n(:,3).^2);

xyz = [ ( -n(:,1).*n(:,3).*neu(:,1) - n(:,2).*neu(:,2) )./cphi + n(:,1).*neu(:,3)   ...
        ( -n(:,2).*n(:,3).*neu(:,1) + n(:,1).*neu(:,2) )./cphi + n(:,2).*neu(:,3)   ...
           cphi.*neu(:,1)                                      + n(:,3).*neu(:,3)   ];

if strcmpi(mode,'xr')
   % output cartesian coordinates instead of differences -> add XR
   if m2 == 1 
     xyz=xyz+repmat(xr,[size(xyz,1) 1]);
   else
     xyz=xyz+xr;
  end
end

if n~=3 & m==3, xyz=xyz';, end

% Optionally output the rotation matrices

if nargout > 1

  %   R(;,:,k) = [ -n(k,1).*n(k,3)./cphi   -n(k,2)./cphi   n(k,1) ;  ...
  %                -n(k,2).*n(k,3)./cphi    n(k,1)./cphi   n(k,2) ; ...
  %                 cphi                    0              n(k,3) ]
  %
  %   Remember, in Matlab, the first index runs fastest 

  l=size(n,1);
  R = [ -n(:,1).*n(:,3)./cphi  -n(:,2).*n(:,3)./cphi  cphi   ...    
        -n(:,2)./cphi           n(:,1)./cphi          zeros(l,1)  ...
         n(:,1)                 n(:,2)                n(:,3) ];     % vector representation of matrix (columns run fastest)
  R=reshape(R',[3 3 l]);                                            % another transpose is needed to correctly align rows of n into columns
        
end


return
