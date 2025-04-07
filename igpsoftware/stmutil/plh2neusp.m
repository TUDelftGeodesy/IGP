function [neu,plh0] = plh2neusp(plh,plh0)
%PLH2NEUSP   Ellipsoidal (Lat,Lon,Hgt) to local coordinates (North,East,Up).
%   [NEU,PLH0]=PLH2NEUSP(PLH) converts the N-by-3 matrix PLH with in the rows 
%   ellipsoidal coordinates Phi, Lambda and h into a N-by-3 matrix NEU
%   with local coordinates (North, East, Up) with respect to a reference
%   point with ellipsoidal coordinates PLH0. The latitude Phi and 
%   longitude Lambda in PLH and PLH0 are in degrees, h in PLH and
%   PLHREF is in meters, and the local coordinates NEU are in meters.
%
%   NEU=PLH2NEUSP(PLH,PLH0) allows the user to specify the reference 
%   point. 
%
%   This function approximates a stereographic double projection (line RD).
%
%   Examples:                                   % plh and plhref in
%       [neu,plh0]=plh2neusp(plh)               % degrees/meter
%     
%   See also NEU2PLHSP, PLH2NEU, XYZ2NEU, XYZ2PLH, XYZ2ZAS, NEU2XYZ, ZAS2XYZ and PLH2XYZ.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016,2020.

%   Created:    30 June 2016 by Hans van der Marel
%   Modified:   28 Sept 2020 by Hans van der Marel
%                - modified PLH2NEU for approximate stereographic projection 
%                - agrees within a few meters with Dutch RD projection in
%                  distances (because of the ellipsoid), coordinate differences 
%                  are larger due to the rotation of RD

% Input argument checking

if nargin < 1
  error('Must be called with at least one argument.');
end
if size(plh,2) ~=3
  error('plh must have three columns or elements')
end

% Reference point latitude (rad), longitue (rad), height (m)

if nargin < 2
   % reference point is in the center, but we keep the original height, ie reference height is zero 
   plh0=[ (min(plh(:,1))+max(plh(:,1)))./2 (min(plh(:,2))+max(plh(:,2)))./2  0];
else
   plh0=plh0(:)';
   if size(plh0,1) == 2 
      plh0=[plh0 0];
   elseif size(plh0,2) ~= 3
      error('plh must have two or three elements')
   end
end

% Ellipsoidal parameters

a=6378137.;         % GRS80(WGS84)
f=1/298.257223563;  % GRS80
e2=2*f-2*f^2;

% Curvatures

N = a ./ sqrt(1 - e2 .* sind(plh0(1)).^2);
M = N * (1 -e2) / ( 1 - e2 .* sind(plh0(1)).^2 );

% convert to NEU

neu = [ (plh(:,1)-plh0(1))*M*pi/180 (plh(:,2)-plh0(2)).*(cosd(plh(:,1))+cosd(plh0(1)))./2*N*pi/180 plh(:,3)-plh0(3) ]; 

end
