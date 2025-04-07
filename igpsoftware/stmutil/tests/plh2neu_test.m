%% Test of plh2neusp and neu2plhsp
%

%% Simulate a regular (lat/lon) grid

lat1=51:.1:52;
lon1=4:.1:5;

[lat,lon] = meshgrid(lat1,lon1);
h=rand(size(lat));

%% Convert lat/lon to neu and compute distance from center point

[neu,plh0]=plh2neusp([lat(:) lon(:) h(:)]);
distneu = sqrt( neu(:,1).^2 + neu(:,2).^2 );

figure
subplot(1,2,1)
plot(lon(:),lat(:),'*')
axis equal
subplot(1,2,2)
plot(neu(:,2)./1000,neu(:,1)./1000,'*')
axis equal

%% Test inverse transformation

plhinv=neu2plhsp(neu,plh0);

max(abs([lat(:) lon(:) h(:)]-plhinv))


%% Compute rd coordinates as reference

rd=etrs2rd([lat(:) lon(:)],0); 
rd0=etrs2rd([plh0(1) plh0(2)],0);

rd(:,1)=rd(:,1)-rd0(1);
rd(:,2)=rd(:,2)-rd0(2);

distrd = sqrt( rd(:,1).^2 + rd(:,2).^2 );

%% Compute distance from stereodp

mstruct.origin(1)=plh0(1);
mstruct.origin(2)=plh0(2);
mstruct.a=6378137.;         % GRS80(WGS84)
mstruct.inv_f=298.257223563;  % GRS80;
mstruct.scale=1;
mstruct.falseeasting=0;
mstruct.falsenorthing=0;


[x,y]=stereodp(lat(:)*pi/180,lon(:)*pi/180,mstruct,'forward');

distsp=sqrt(x.^2+y.^2);
figure
plot(distneu-distsp,'linewidth',1)
hold on
plot(neu(:,1)-y)
plot(neu(:,2)-x)
legend('dist','x','y')


%% Compute distance from sperical coordinates

a=6378137.;         % GRS80(WGS84)
f=1/298.257223563;  % GRS80
e2=2*f-2*f^2;
b=a*(1-f);

phi=mean(lat);
Re=sqrt( ( (a.^2*cosd(phi)).^2+(b.^2*sind(phi)).^2 ) / ( (a.*cosd(phi)).^2+(b.*sind(phi)).^2 ) );

dists = distance(lat(:),lon(:),plh0(1),plh0(2),Re);

%% Evaluate the differences

figure
plot(distneu-distrd,'linewidth',1)
%hold on
%plot(distneu-dists)
%plot(dists-distrd)
%legend('neu-rd','neu-s','s-rd')

d=reshape(distneu-distrd,size(lon));

figure
surf(lat1,lon1,d,d)
colorbar

figure
contour(lon1,lat1,d')
colorbar

%% Local functions

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
%   (c) Hans van der Marel, Delft University of Technology, 2016.

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

function plh = neu2plhsp(neu,plh0)
%NEU2PLHSP   Local coordinates (North,East,Up) to Ellipsoidal (Lat,Lon,Hgt).
%   PLH=NEU2PLH2SP(PLH) converts the N-by-3 matrix NEU with in the rows
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
%   (c) Hans van der Marel, Delft University of Technology, 2016.

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
   