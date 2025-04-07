function xyzetrs=rdnap2etrs(rdnap,varargin)
%RDNAP2ETRS Convert RDNAP coordinates into cartesian ETRS89 coordinates.
%  Convert RDNAP coordinates into cartesian ETRS89 coordinates.
%  Syntax
%             xyzetrs=rdnap2etrs(rdnap)
%             xyzetrs=rdnap2etrs(rdnap,options)
%  Input
%    rdnap    n-by-3 array with RDNAP coordiantes [x_rd,y_rd,nap]
%  Options
%    id       system [MOSTRECENT|RDNAPTRANS2008|RDNAPTRANS2004], default most 
%             recent
%    outcrd   output coordinates [XYZ|PLH], default cartesian XYZ coordinates
%    extrd    RD correction grid extension [NaN|zero], default NaN
%  Output
%    xyzetrs  n-by-3 array with cartesian coordinates [X,Y,Z] or [lat,lon,hell]
%             depending on the optional argument outcrd
% 
%  See also etrs2rdnap, etrs2nap and nap2etrs.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013

% Created:   6 Dec 2004 by Hans van der Marel, TUD
% Modified: 22 Apr 2010 by Sierd de Vries
%              - Added NAP correction for RDNAPTRANS2008
%           23 May 2011 by Piers Titus van der Torren
%              - Vectorized the function
%            6 Jun 2013 by Hans van der Marel 
%              - integrated modifications of 22 Apr 2010 and 23 May 2011 with
%                some modifications
%              - Select between RDNAPTRANS2004 and RDNAPTRANS2008 
%              - Call nap2etrs to do height transformation instead of hardwired
%                transformation with correction
%              - Added option to output latitude, longitude and height
%                instead of cartesian coordinates, and extend RD coordinates
%           28 Mar 2018 by Hans van der Marel
%              - Added check with 'ZERO' option to prevent NaN due to invalid height

MEAN_GEOID_HEIGHT_BESSEL = 0;

%  Check the input arguments

if (nargin < 1)
  error('must have at least 1 input argument')
end

% process the options

id='MOSTRECENT';
outcrd='XYZ';
extrd='NAN';

for i=1:length(varargin)
  switch upper(varargin{i})
     case {'RDNAPTRANS2008','2008','RDNAPTRANS2004','2004','MOSTRECENT','CURRENT'}
       id=varargin{i};
     case {'PLH','GEO','XYZ'}
       outcrd=varargin{i};
     case {'NAN','ZERO'}
       extrd=varargin{i};
     otherwise
       error('invalid option')
  end
end

% Check size of input array (for single coordinate input allow column vector)

if all(size(rdnap) == [3,1])
  rdnap = rdnap.';
  is3x1 = true;
else
  is3x1 = false;
end
assert(size(rdnap,2)==3,'input array must have shape n-by-3')

%  Upon first use, load the correction tables from the x2c.grd and
%  y2c.grd files

global RDGRIDX RDGRIDY
if isempty(RDGRIDX) || isempty(RDGRIDY)
  RDGRIDX=readgrd('x2c.grd');
  RDGRIDY=readgrd('y2c.grd');
end

% Input parameters

xrd=rdnap(:,1);
yrd=rdnap(:,2);
nap=rdnap(:,3);

% Calculate appoximated value of ellipsoidal Bessel height
% The error made by using a constant for de Bessel geoid height is 
% circa 1 meter in the ellipsoidal height (for the NLGEO2004 geoid). 
% This intoduces an error in the phi, lambda position too, this error is 
% nevertheless certainly smaller than 0.0001 m.
% 
h_bessel = nap + MEAN_GEOID_HEIGHT_BESSEL;
h_bessel(isnan(nap))=MEAN_GEOID_HEIGHT_BESSEL;

% Compute corrections for the RD coordinates

dx=grdint(RDGRIDX,xrd,yrd);
dy=grdint(RDGRIDY,xrd,yrd);

if strcmpi(extrd,'ZERO')
   dx(isnan(dx))=0;
   dy(isnan(dy))=0;
end

xpseudord=xrd+dx;
ypseudord=yrd+dy;

[phi_bessel,lam_bessel]=stereodp(xpseudord, ypseudord,nlbessel,'inverse');

% Coordinates using approximate height

plhbessel=[phi_bessel,lam_bessel,h_bessel];

xyzbessel=plh2xyz(plhbessel,'BESSEL');
xyzetrs=simtrans(xyzbessel,bessel2etrs(id));

% Take the correct height into account

plhetrs=xyz2plh(xyzetrs,'GRS80');
%plhetrs(:,3)=nap2etrs(plhetrs(:,1)*180/pi,plhetrs(:,2)*180/pi,nap,id);
hell=nap2etrs(plhetrs(:,1)*180/pi,plhetrs(:,2)*180/pi,nap,id);
if strcmpi(extrd,'ZERO')
   validhgt=find(~isnan(hell));
   plhetrs(validhgt,3)=hell(validhgt);
else
   plhetrs(:,3)=hell;
end

switch upper(outcrd)
   case {'XYZ'}
      xyzetrs=plh2xyz(plhetrs,'GRS80');
   case {'PLH','GEO'} 
      xyzetrs=plhetrs;
end

if is3x1
   xyzetrs=xyzetrs';
end

return;
