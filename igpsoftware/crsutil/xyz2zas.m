function [z,a,s] = xyz2zas(xyz,plh,mode)
%XYZ2ZAS   Compute zenith angle, azimuth and distance from cartesian coordinates.
%   [Z,A,S]=XYZ2ZAS(DXYZ,PLH) computes the M-by-1 vectors Z, A, and S with 
%   the zenith angle, azimuth and distance from the M-by-3 matrix with cartesian
%   coordinate differences (dX, dY and dZ) in the ECEF system. PLH contains the 
%   latitude and longitude (height is not used), in radians, corresponding to 
%   the zenith direction. PLH is either a M-by-2 or M-by-3 matrix, or a 1-by-2 
%   or 1-by-3 vector. In the latter case the same center point is used for all 
%   coordinate triplets.
%
%   This function has an optional third argument MODE. Default for MODE is 'p'
%   ('plh'). Other MODE's are 'n' ('normal'), 'r' and 'xr':
%
%   [Z,A,S]=XYZ2ZAS(DXYZ,N,'n') uses the normal vector N instead of PLH. N is a 
%   N-by-3 matrix or 1-by-3 vector vector pointing in the up direction (N does 
%   not have to be a unit vector).
%
%   [Z,A,S]=XYZ2ZAS(DXYZ,XR,'r') uses the cartesian coordinates XR of the center 
%   point instead of geographic coordinates PLH. XR is a N-by-3 matrix or 
%   1-by-3 vector. 
%
%   [Z,A,S]=XYZ2ZAS(XYZ,XR,'xr') uses cartesian coordinates XYZ (instead of 
%   coordinate differences DXYZ) and the center point XR as with the previous 
%   option. 
%
%   Examples:
%       [z,a,s]=xyz2zas(dxyz,[52*pi/180 4*pi/180])
%       [z,a,s]=xyz2zas(dxyz,plh)
%       [z,a,s]=xyz2zas(dxyz,n,'n')
%       [z,a,s]=xyz2zas(dxyz,xr,'r')
%       [z,a,s]=xyz2zas(xyz,xr,'xr')
%
%   Please note that the following two examples do not give the same results
%   because the up-direction is different (geodetic versus astronomic)
%
%       xr=plh2xyz(plh) 
%       [z,a,s]=xyz2zas(dxyz,xr,'r')  (same result as [z,a,s]=xyz2zas(dxyz,plh))
%       [z,a,s]=xyz2zas(dxyz,xr,'n')
%
%   See also XYZ2NEU, XYZ2PLH, NEU2XYZ, ZAS2XYZ and PLH2XYZ.
%
%   (c) Hans van der Marel, Delft University of Technology, 2007, 2013.

%   Created:    24 May 2007 by Hans van der Marel
%   Modified:   14 Jun 2013 by Hans van der Marel
%               - updated description

% Input argument checking

if nargin < 2 | nargin > 3
  error('Must be called by two or three arguments.');
end
if nargin == 2
  mode='plh';
end

[m,n]=size(xyz);
if n~=3 & m==3
  xyz=xyz';
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
   otherwise
      error('unknown mode')
end

% Do the actual computation 

%
% neu = [ ( -n(:,1).*xyz(:,1) - n(:,2).*xyz(:,2) ) .* n(:,3) ./ cphi + cphi.*xyz(:,3)     ...
%         ( -n(:,2).*xyz(:,1) + n(:,1).*xyz(:,2) ) ./  cphi                              ...
%            n(:,1).*xyz(:,1) + n(:,2).*xyz(:,2) + n(:,3).*xyz(:,3) ]
%
% cphi= sqrt(1-n(:,3).^2);
% ip=n(:,1).*xyz(:,1) + n(:,2).*xyz(:,2) + n(:,3).*xyz(:,3);
%
% neu = [ (  ip .* -n(:,3)  + xyz(:,3) ) ./ cphi               ...
%         ( -n(:,2).*xyz(:,1) + n(:,1).*xyz(:,2) ) ./  cphi    ...
%            ip                                                  ];
%
% s = sqrt( xyz(:,1).*xyz(:,1) + xyz(:,2).*xyz(:,2) + xyz(:,3).*xyz(:,3) );
% z = acos(neu(:,3)./s);
% a = atan2(neu(:,2),neu(:,1));

ip=n(:,1).*xyz(:,1) + n(:,2).*xyz(:,2) + n(:,3).*xyz(:,3);

s=  sqrt( xyz(:,1).*xyz(:,1) + xyz(:,2).*xyz(:,2) + xyz(:,3).*xyz(:,3) );
z = acos(ip./s);
a = atan2(  -n(:,2).*xyz(:,1) + n(:,1).*xyz(:,2) ,   ip .* -n(:,3)  + xyz(:,3) );

return
