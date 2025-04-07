function plh = neu2plhsp(neu,plh0)
%NEU2PLHSP   Local coordinates (North,East,Up) to Ellipsoidal (Lat,Lon,Hgt).
%   PLH=NEU2PLH2SP(NEU,PLH0) converts the N-by-3 matrix NEU with in the rows
%   local coordinates (North, East, Up), with respect to a reference point 
%   with ellipsoidal coordinates PLH0, into a N-by-3 matrix PLH with 
%   ellipsoidal coordinates Phi, Lambda and h. The latitude Phi and 
%   longitude Lambda in PLH and PLH0 are in degrees, h in PLH and
%   PLHREF is in meters, and the local coordinates NEU are in meters.
%
%   This function approximates a stereographic double projection (like RD).
%
%   Examples:                                   % plh and plhref in
%       plh=neu2plhsp(neu,plh0)                 % degrees/meter
%     
%   See also PLH2NEUSP, PLH2NEU, XYZ2NEU, XYZ2PLH, XYZ2ZAS, NEU2XYZ, ZAS2XYZ and PLH2XYZ.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016,2020.

%   Created:    28 Sept 2020 by Hans van der Marel
%   Modified:   

% Input argument checking

if nargin ~=2
  error('Must be called with two input arguments.');
end
if size(neu,2) ~=3
  error('neu must have three columns or elements')
end

% Reference point latitude (rad), longitue (rad), height (m)

plh0=plh0(:)';
if size(plh0,1) == 2 
   plh0=[plh0 0];
elseif size(plh0,2) ~= 3
  error('plh must have two or three elements')
end

% Ellipsoidal parameters

a=6378137.;         % GRS80(WGS84)
f=1/298.257223563;  % GRS80
e2=2*f-2*f^2;

% Curvatures

N = a ./ sqrt(1 - e2 .* sind(plh0(1)).^2);
M = N * (1 -e2) / ( 1 - e2 .* sind(plh0(1)).^2 );

% convert to NEU to PLH

lat = plh0(1) + neu(:,1)/(M*pi/180); 
lon = plh0(2) + neu(:,2)/(N*pi/180) ./ ((cosd(lat)+cosd(plh0(1)))./2); 

plh= [ lat lon +plh0(3)+neu(:,3) ]; 
    
end
   