function hnap=etrs2nap(phi,lam,hell,id)
%ETRS2NAP Convert ellipsoidal ETRS89 height into NAP height
%  Convert ellipsoidal ETRS89 height into NAP height
%  Syntax
%           hnap=etrs2nap(phi,lam,hell)
%           hnap=etrs2nap(phi,lam,hell,id)
%  Input
%    phi    lattitude (degrees)
%    lam    longitude (degrees)
%    hell   height above ETRS89 ellipsoid (m)
%    id     optional id [RDNAPTRANS2008|RDNAPTRANS2004], default most recent  
%  Output
%    hnap   NAP height
% 
%  See also rdnap2etrs, etrs2rdnap and nap2etrs.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013
    
%  Created:   6 Dec 2004 by Hans van der Marel, TUD
%  Modified   6 Jun 2013 by Hans van der Marel
%               - Select between RDNAPTRANS2004 and RDNAPTRANS2008 

% Check the input arguments

if (nargin < 3)
   error('must have at least 3 input arguments')
elseif (nargin < 4) 
   id='MOSTRECENT';
end

switch upper(id)
   case {'RDNAPTRANS2008','2008','MOSTRECENT','CURRENT'}
     dhnap=-.0088;
   case {'RDNAPTRANS2004','2004'}
     dhnap=0;
   otherwise
     error(['Unknown system ' id ])
end

% Upon first use, load the NLGEOID from the nlgeo04.grd file

global NLGEOID
if isempty(NLGEOID)
  NLGEOID=readgrd('nlgeo04.grd');
end

% Interpolate the geoid height N

N=grdint(NLGEOID,lam,phi);

% Compute the NAP height

hnap=hell-N-dhnap;

return;

