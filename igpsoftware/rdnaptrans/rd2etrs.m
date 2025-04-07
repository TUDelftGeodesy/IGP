function latlon=rd2etrs(rd,hmsl)
%RD2ETRS  Convert RD coordinates into latitude and longitude.
%  LATLON=RD2ETRS(RD) convert matrix RD with Dutch RD coordinates into
%  output LATLON with ETRS89 latitude and longitude in degrees. LATLON is 
%  a matrix of the same size as RD with in the first column the ETRS89 
%  latitude in degrees, and in its second column the longitude in degrees. 
%  RD contains the grid coordinates in meters. The height is assumed to be
%  zero in the Dutch NAP system.
%
%  LATLON=RD2ETRS(RD,HMSL) does the same, but uses HMSL for the height
%  in the Dutch NAP system. This height does not have to be very accurate,
%  and has only an effect at the mm level. HMSL can be scalar or an array.
%
%  Example:
%  
%    latlon=etrs2rd([ 193371.413   308271.466 ]);
% 
%  See also etrs2rd, etrs2rdnap, rdnap2etrs, etrs2nap and nap2etrs.
% 
%  (c) Hans van der Marel, Delft University of Technology, 2004-2019

% Created:  30 Apr 2019 by Hans van der Marel, TUD
% Modified:  1 Oct 2023 by Hans van der Marel, TUD
%            - Set RD correction grid extension to ZERO to avoid NaN's

% Check input arguments

if nargin < 1 || nargin > 2, error('Function needs one or two input arguments.'); end
if nargin < 2 
   hmsl=0;
end

if isscalar(hmsl)
  hmsl=repmat(hmsl,[size(rd,1),1]);
end    

% Compute latitude and longitude

xyz=rdnap2etrs([ rd hmsl],'ZERO');
plh=xyz2plh(xyz);

latlon=plh(:,1:2)*180/pi;

end