function rdnap=etrs2rdnap(xyzetrs,varargin)
%ETRS2RDNAP Convert cartesian ETRS89 coordinates into RDNAP.
%  Convert cartesian ETRS89 coordinates in RDNAP coordinates.
%  Syntax
%             rdnap=etrs2rdnap(xyzetrs)
%             rdnap=etrs2rdnap(xyzetrs,options)
%  Input
%    xyzetrs  n-by-3 array with cartesian coordinates [X,Y,Z] or [lat,lon,hell] 
%             depending on the optional argument incrd (units meters/radians) 
%  Options
%    id       system [MOSTRECENT|RDNAPTRANS2008|RDNAPTRANS2004], default most 
%             recent
%    incrd    input coordinates [XYZ|PLH], default cartesian XYZ coordinates
%    extrd    RD correction grid extension [NaN|zero], default NaN
%  Output
%    rdnap    n-by-3 array with RDNAP coordiantes [x_rd,y_rd,nap]
% 
%  Examples:
%    rdnap=etrs2rdnap([3924697.7800 301125.1300 5001905.2900]);
%    rdnap=etrs2rdnap([51.8237*pi/180    5.7931*pi/180   52.1970],'PLH');
%
%  See also rdnap2etrs, etrs2nap and nap2etrs.
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
%              - Call etrs2nap to do height transformation instead of hardwired
%                transformation with correction
%              - Added option to input latitude, longitude and height
%                instead of cartesian coordinates, and extend RD coordinates

%  Check the input arguments

if (nargin < 1)
  error('must have at least 1 input argument')
end

% process the options

id='MOSTRECENT';
incrd='XYZ';
extrd='NAN';

for i=1:length(varargin)
  switch upper(varargin{i})
     case {'RDNAPTRANS2008','2008','RDNAPTRANS2004','2004','MOSTRECENT','CURRENT'}
       id=varargin{i};
     case {'PLH','GEO','XYZ'}
       incrd=varargin{i};
     case {'NAN','ZERO'}
       extrd=varargin{i};
     otherwise
       error('invalid option')
  end
end


% Check size of input array (for single coordinate input allow column vector)

if all(size(xyzetrs) == [3,1])
  xyzetrs = xyzetrs.';
  is3x1 = true;
else
  is3x1 = false;
end
assert(size(xyzetrs,2)==3,'input array must have shape n-by-3')

%  Upon first use, load the correction tables from the x2c.grd and
%  y2c.grd files

global RDGRIDX RDGRIDY
if isempty(RDGRIDX) || isempty(RDGRIDY)
  RDGRIDX=readgrd('x2c.grd');
  RDGRIDY=readgrd('y2c.grd');
end

% Compute phi, lambda and height above the Bessel ellipsoid

switch upper(incrd)
   case {'PLH','GEO'} 
      invalidhgt=find(isnan(xyzetrs(:,3)));
      xyzetrs(invalidhgt,3)=0;
      xyzetrs=plh2xyz(xyzetrs,'GRS80');
   otherwise
      invalidhgt=[];
end
xyzbessel=simtrans(xyzetrs,etrs2bessel(id));
plhbessel=xyz2plh(xyzbessel,'BESSEL');   

% Compute pseudo-rd coordinates and corrections for RD coordinates

[xpseudord,ypseudord]=stereodp(plhbessel(:,1),plhbessel(:,2),nlbessel,'forward');

dx=grdint(RDGRIDX,xpseudord,ypseudord);
dy=grdint(RDGRIDY,xpseudord,ypseudord);

if strcmpi(extrd,'ZERO')
   dx(isnan(dx))=0;
   dy(isnan(dy))=0;
end

xrd=xpseudord-dx;
yrd=ypseudord-dy;

% Compute NAP height (using lattitude and longitude in ETRS)

plhetrs=xyz2plh(xyzetrs,'GRS80') ;
nap=etrs2nap(plhetrs(:,1)*180/pi,plhetrs(:,2)*180/pi,plhetrs(:,3),id);

% Set output variable

nap(invalidhgt)=nan;
rdnap=[xrd,yrd,nap];

if is3x1
   rdnap=rdnap';
end

return;
