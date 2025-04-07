function xyz = plh2xyz(plh,ellips)
%PLH2XYZ Ellipsoidal coordinates to Cartesian Coordinates 
%   XYZ=PLH2XYZ(PLH) converts the Nx3 matrix PLH with in the rows ellipsoidal 
%   coordinates Phi, Lambda and h into a Nx3 matrix XYZ with cartesian 
%   coordinates X, Y and Z. Phi and Lambda are in radians, h is in meters.
%
%   XYZ=PLH2XYZ(PLH,ELLIPS) allows to specify the ellipsoid. ELLIPS is a text
%   is a text string with the name of the ellipsoid or a vector with the 
%   semi-major axis a and flattening 1/f. Default for ellips is 'WGS-84'.
%
%   See also XYZ2PLH, PLH2STR, STR2PLH, XYZ2NEU, NEU2XYZ and INQELL.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995,2013

%   Created:     7 May 1995 by Hans van der Marel
%   Modified:   14 Jun 2013 by Hans van der Marel
%                - updated description

if nargin<2, ellips='WGS-84';, end
if isstr(ellips)
  [a,f] = inqell(ellips);
else
  a=ellips(1);
  f=1/ellips(2);
end
% excentricity e (squared) 
e2 = 2*f - f^2;

[m,n]=size(plh);
if n~=3 & m==3, plh=plh';, end

N = a ./ sqrt(1 - e2 .* sin(plh(:,1)).^2);
xyz = [ (N+plh(:,3)).*cos(plh(:,1)).*cos(plh(:,2)) ...
        (N+plh(:,3)).*cos(plh(:,1)).*sin(plh(:,2)) ...
        (N-e2.*N+plh(:,3)).*sin(plh(:,1))          ];

if n~=3 & m==3, xyz=xyz';, end
