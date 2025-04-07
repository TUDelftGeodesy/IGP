function xyz = zas2xyz(z,a,s,plh,mode)
%ZAS2XYZ Zenith angle, azimuth and distance to cartesian coordinates (dX,dY,dZ).
%   DXYZ=ZAS2XYZ(Z,A,S,PLH) computes the M-by-3 matrix DXYZ with cartesian 
%   coordinate differences (dX, dY and dZ), in a ECEF system, from the M-by-1
%   vectors Z, A and S with the zenith angle, geodetic azimuth and distance. 
%   PLH contains the latitude and longitude (height is not used), in radians, 
%   corresponding to the zenith direction. PLH is either a M-by-2 or M-by-3 
%   matrix, or a 1-by-2 or 1-by-3 vector. In the latter case the same center 
%   point is used for all coordinate triplets.
%
%   This function has an optional fifth argument MODE. Default for MODE is 'p'
%   ('plh'). Other MODE's are 'n' ('normal'), 'r' and 'xr':
%
%   DXYZ=ZAS2XYZ(Z,A,S,N,'n') uses the normal vector N instead of geographic 
%   coordinates PLH of the center point. N is a N-by-3 matrix or 1-by-3 vector 
%   vector pointing in the up direction (N does not have to be a unit vector).
%
%   DXYZ=ZASXYZ(Z,A,S,XR,'r') uses the cartesian coordinates XR of the center 
%   point instead of geographic coordinates PLH. XR is a N-by-3 matrix or 
%   1-by-3 vector. 
%
%   XYZ=ZAS2XYZ(A,A,S,XR,'xr') outputs cartesian coordinates XYZ (instead of 
%   coordinate differences DXYZ), and uses the center point XR as with the 
%   previous option. XYZ is computed inside the function from XR.
%
%   Examples:
%       dxyz=zas2xyz(z,a,s,[52*pi/180 4*pi/180])
%       dxyz=zas2xyz(z,a,s,plh)
%       dxyz=zas2xyz(z,a,s,n,'n')
%       dxyz=zas2xyz(z,a,s,xr,'r')
%       xyz=zas2xyz(z,a,s,xr,'xr')
%
%   Please note that the following two examples do not give the same results
%   because the up-direction is different (geodetic versus astronomic)
%
%       xr=plh2xyz(plh) 
%       dxyz=zas2xyz(z,a,s,xr,'r')    (same result as dxyz=zas2xyz(z,a,s,plh))
%       dxyz=zas2xyz(z,a,s,xr,'n')
%
%   See also NEU2XYZ, XYZ2ZAS, XYZ2NEU, XYZ2PLH and PLH2XYZ.
%
%   (c) Hans van der Marel, Delft University of Technology, 2013.

%   Created:    14 Jun 2013 by Hans van der Marel
%   Modified:   

% Input argument checking

if nargin < 4 | nargin > 5
  error('Must be called by four or five arguments.');
end
if nargin == 4
  mode='plh';
end

[m,n]=size(z);
if n~=1 & m==1
  z=z';a=a';s=s';
end
m1=size(z,1);
if size(a,1) ~= m1 | size(s,1) ~= m1
   error(['Z, A and S must have the same length']);
end
if size(z,2) ~= 1 | size(a,2) ~= 1 | size(s,2) ~= 1
   error(['Z, A and S must be vectors']);
end
[m2,n2]=size(plh);
if n2~=3 & n2~=2 & ( m2==3 | m2==2 )
  plh=plh'; 
  [m2,n2]=size(plh);
end
if m1~=m2 & m2~=1
   error(['PLH, N or XR must be same length as Z, A and S, or, a vector']);
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
   otherwise
      error('unknown mode')
end

% Do the actual computation 

% xyz = [  -n(:,1).*n(:,3)./cphi.*neu(:,1)  - n(:,2)./cphi.*neu(:,2) + n(:,1).*neu(:,3)   ...
%          -n(:,2).* n(:,3)./cphi.*neu(:,1) + n(:,1)./cphi.*neu(:,2) + n(:,2).*neu(:,3)   ...
%           cphi.*neu(:,1)                                           + n(:,3).*neu(:,3)    ]

cphi= sqrt(1-n(:,3).^2);

neu= [ s.*cos(a).*sin(z)  s.*sin(a).*sin(z)  s.*cos(z) ];
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

if n~=1 & m==1, xyz=xyz';, end

return
