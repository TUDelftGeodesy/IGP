function [neu,R] = xyz2neu(xyz,plh,mode)
%XYZ2NEU   Cartesian coordinates (dX,dY,dZ) to local coordinates (North,East,Up).
%   NEU=XYZ2NEU(DXYZ,PLH) converts the M-by-3 matrix DXYZ with rows of cartesian
%   coordinate differences (dX, dY and dZ), ECEF system, into the M-by-3 matrix
%   NEU with local coordinates (North, East, Up). PLH contains the geographic 
%   coordinates, i.e. latitude and longitude (height is not used) in radians,
%   of the center of the local system. PLH is either a M-by-2 or M-by-3 matrix, 
%   or a 1-by-2 or 1-by-3 vector. In case PLH is a vector the same center point 
%   is used for all coordinate triplets.
%
%   This function has an optional third argument MODE. Default for MODE is 'p'
%   ('plh'). Other MODE's are 'n' ('normal'), 'r' and 'xr':
%
%   NEU=XYZ2NEU(DXYZ,N,'n') uses the normal vector N instead of geographic 
%   coordinates PLH of the center point. N is a M-by-3 matrix or 1-by-3 vector 
%   vector pointing in the up direction (N does not have to be a unit vector).
%
%   NEU=XYZ2NEU(DXYZ,XR,'r') uses the cartesian coordinates XR of the center 
%   point instead of geographic coordinates PLH. XR is a M-by-3 matrix or 
%   1-by-3 vector. 
%
%   NEU=XYZ2NEU(XYZ,XR,'xr') uses cartesian coordinates XYZ (instead of 
%   coordinate differences DXYZ) and the center point XR as with the previous 
%   option. DXYZ is computed inside the function from XR.
%
%   [NEU,R]=XYZ2NEU(...) also outputs the rotation matrix R, with NEU'=R*DXYZ'.
%   R is a 3-by-3 matrix in case PLH, N or XR is a vector. In case PLH, N,
%   or XR is a matrix, R is 3-by-3-by-M. The rotation matrix can be used for 
%   instance to convert the covariance matrix (when available and needed).
%
%   Examples:
%       neu=xyz2neu(dxyz,[52*pi/180 4*pi/180])
%       neu=xyz2neu(dxyz,plh)
%       neu=xyz2neu(dxyz,n,'n')
%       neu=xyz2neu(dxyz,xr,'r')
%       neu=xyz2neu(xyz,xr,'xr')
%       [neu,R]=xyz2neu(xyz,xr,'xr')
%
%   Please note that the following two examples do not give the same results
%   because the up-direction is different (geodetic versus astronomic)
%
%       xr=plh2xyz(plh) 
%       neu=xyz2neu(dxyz,xr,'r')    (same result as neu=xyz2neu(dxyz,plh))
%       neu=xyz2neu(dxyz,xr,'n')
%
%   See also XYZ2PLH, XYZ2ZAS, NEU2XYZ, ZAS2XYZ and PLH2XYZ.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995,2007,2013.

%   Created:     7 May 1995 by Hans van der Marel
%   Modified:   24 May 2007 by Hans van der Marel
%                - faster implementation reducing use of trigonometry 
%                - option to use normal vector N instead of phi/lambda
%                - updated description
%               14 Jun 2013 by Hans van der Marel
%                - minor changes to description
%                - added the options to use XR
%               12 March 2014 by Hans van der Marel
%                - added the option to output the rotation matrix
%                - minor corrections to description
%               12 April 2020 by Hans van der Marel
%                - more robust test on transpose (for Matlab 2020?)

% Input argument checking

if nargin < 2 | nargin > 3
  error('Must be called by two or three arguments.');
end
if nargin == 2
  mode='plh';
end

[m,n]=size(xyz);
if n~=3 & m==3
  dotranspose=true;
  xyz=xyz';
else
  dotranspose=false;
end
[m1,n1]=size(xyz);
[m2,n2]=size(plh);
if n2~=3 & n2~=2 & ( m2==3 | m2==2 )
  plh=plh'; 
  [m2,n2]=size(plh);
end
if m1~=m2 & m2~=1
   error(['PLH, N or XR must be same length as XYZ or a vector']);
end

% Sort out the various modes

switch lower(mode)
   case {'normal','n'}
      % plh contains normal vector -> normalize to unit vector
      n=plh ./ sqrt( plh(:,1).*plh(:,1) + plh(:,2).*plh(:,2) + plh(:,3).*plh(:,3) );
   case {'plh','p'}
      % plh contains latutude and longitude -> compute unit normal vector
      n= [ cos(plh(:,1)).*cos(plh(:,2)) cos(plh(:,1)).*sin(plh(:,2))  sin(plh(:,1)) ];
   case {'r'}
      % plh contains XR -> compute latitude and longitude, then unit normal vector
      plh=xyz2plh(plh);
      n= [ cos(plh(:,1)).*cos(plh(:,2)) cos(plh(:,1)).*sin(plh(:,2))  sin(plh(:,1)) ];
   case {'xr'}
      % xyz contains cartesian coordinates -> compute coordinate differences with XR
      if m2 == 1 
         xyz=xyz-repmat(plh,[size(xyz,1) 1]);
      else
         xyz=xyz-plh;
      end
      % plh contains XR -> compute latitude and longitude, then unit normal vector
      plh=xyz2plh(plh);
      n= [ cos(plh(:,1)).*cos(plh(:,2)) cos(plh(:,1)).*sin(plh(:,2))  sin(plh(:,1)) ];
   case {'legacy'}
      % Old way of doing things (using less efficient and flexible trigonometry)
      neu = [ -sin(plh(:,1)).*cos(plh(:,2)).*xyz(:,1) - sin(plh(:,1)).*sin(plh(:,2)).*xyz(:,2) + cos(plh(:,1)).*xyz(:,3) ...
                             -sin(plh(:,2)).*xyz(:,1) +                cos(plh(:,2)).*xyz(:,2)                           ...
               cos(plh(:,1)).*cos(plh(:,2)).*xyz(:,1) + cos(plh(:,1)).*sin(plh(:,2)).*xyz(:,2) + sin(plh(:,1)).*xyz(:,3) ]
      if n~=3 & m==3, neu=neu';, end
      return
   otherwise
      error('unknown mode')
end

% Do the actual computation 

% neu = [ ( -n(:,1).*xyz(:,1) - n(:,2).*xyz(:,2) ) .* n(:,3) ./ cphi + cphi.*xyz(:,3)     ...
%         ( -n(:,2).*xyz(:,1) + n(:,1).*xyz(:,2) ) ./  cphi                              ...
%            n(:,1).*xyz(:,1) + n(:,2).*xyz(:,2) + n(:,3).*xyz(:,3) ]

cphi= sqrt(1-n(:,3).^2);
ip=n(:,1).*xyz(:,1) + n(:,2).*xyz(:,2) + n(:,3).*xyz(:,3);

neu = [ (  ip .* -n(:,3)  + xyz(:,3) ) ./ cphi               ...
        ( -n(:,2).*xyz(:,1) + n(:,1).*xyz(:,2) ) ./  cphi    ...
           ip                                                  ];

if dotranspose, neu=neu';, end

% Optionally output the rotation matrices

if nargout > 1

  %   R(;,:,k) = [ -n(k,1).*n(k,3)./cphi  -n(k,2).*n(k,3)./ cphi   cphi   ;  ...
  %                -n(k,2)./cphi           n(k,1)./cphi            0      ; ...
  %                 n(k,1)                 n(k,2)                  n(k,3) ]
  %
  %   Remember, in Matlab, the first index runs fastest 

  l=size(n,1);
  R = [ -n(:,1).*n(:,3)./cphi   -n(:,2)./cphi  n(:,1) ...   
        -n(:,2).*n(:,3)./cphi    n(:,1)./cphi  n(:,2) ...
        cphi                     zeros(l,1)    n(:,3) ];     % vector representation of matrix (columns run fastest)
  R=reshape(R',[3 3 l]);                                     % another transpose is needed to correctly align rows of n into columns
        
end

return
