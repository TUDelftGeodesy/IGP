function neu = plh2neu(plh,plhref,varargin)
%PLH2NEU   Ellipsoidal (Lat,Lon,Hgt) to local coordinates (North,East,Up).
%   NEU=PLH2NEU(PLH,PLHREF) converts the N-by-3 matrix PLH with in the rows 
%   ellipsoidal coordinates Phi, Lambda and h into a N-by-3 matrix NEU
%   with local coordinates (North, East, Up) with respect to a reference
%   point with ellipsoidal coordinates PLHREF. The latitude Phi and 
%   longitude Lambda in PLH and PLHREF are in radians, h in PLH and
%   PLHREF is in meters. The local coordinates NEU are in meters.
%
%   NEU=PLH2NEU(PLH,PLHREF,OPTIONS,...) allows to specify options. Valid
%   options are
%
%      'deg','degrees'   lat and lon input in degrees, height in meters
%      'rad','radians'   lat and lon input in radians, height in meters (default)
%      'local'           output in rectangular local North, East, Up system
%      'ell','ellips'    output in curved ellipsoidal North, East, Up system 
%
%   Examples:                                            % plh and plhref in
%       neu=plh2neu(plh,plhref)                          % radians
%       neu=plh2neu([52+1/60 4 40],[52 4 0],'deg')       % degrees
%       neu=plh2neu([52+1/60 4 40],[52 4 0],'deg','ell') % using alternative method
%     
%   See also XYZ2NEU, XYZ2PLH, XYZ2ZAS, NEU2XYZ, ZAS2XYZ and PLH2XYZ.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016.

%   Created:    30 June 2016 by Hans van der Marel
%   Modified:   


% Input argument checking

if nargin < 2
  error('Must be called with at least two arguments.');
end

if size(plh,2) ~=3
  error('plh must have three columns or elements')
end

plhref=plhref(:)';
if size(plhref,1) ~= 1 || size(plhref,2) ~= 3 
  error('plh must have three elements')
end

% Option processing

method=1;

for k=1:numel(varargin)
   switch lower(varargin{k})
      case {'rad','radians'}
         %disp('Input is in radians.')
      case {'deg','degrees'}
         % convert input to radians
         plh(:,1:2)=plh(:,1:2).*pi/180;
         plhref(:,1:2)=plhref(:,1:2).*pi/180;         
      case {'local','cartesian'}
         method=1;          
      case {'ell','ellips'}
         method=2;
      otherwise
         disp([ 'Unknown option ' varargin{k} ])
   end 
end

% convert to NEU

n=size(plh,1);
if method == 1
   xyz=plh2xyz(plh);
   xyzref=plh2xyz(plhref);
   dxyz=xyz-repmat(xyzref,[ n 1]);
   neu=xyz2neu(dxyz,plhref);
elseif method == 2 
   a=6378137.;         % GRS80(WGS84)
   f=1/298.257223563;  % GRS80
   e2=2*f-2*f^2;
   N_curvature=a/sqrt(1-e2*sin(plhref(1,1))^2);
   M_curvature=a*(1-e2)/(1-e2*sin(plhref(1,1))^2)^(3/2);
   dplh=plh-repmat(plhref,[ n 1]);
   neu(:,1)=(M_curvature+plhref(1,3))*dplh(:,1);
   neu(:,2)=(N_curvature+plhref(1,3))*cos(plhref(1,1))*dplh(:,2);
   neu(:,3)=dplh(:,3);
end

end
    
