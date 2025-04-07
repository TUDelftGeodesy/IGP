function crdOut = crstrans(crdIn,crsIn,crsOut,method)
%CRSTRANS    Coordinate reference system transformation.
%   CRDOUT=CRSTRANS(CRDIN,CRSIN,CRSOUT) converts the N-by-2or3 matrix
%   CRDIN with in the rows X,Y or Lat,Lon coordinates and optionally
%   a height H into a N-by-3 matrix CRDOUT. Values in CRDIN should 
%   be in meters or degrees. Values in CRDOUT are also in meters and
%   degrees. CRSIN and CRSOUT are the input and output coordinate
%   reference system, respectively. In principle any geographic or
%   map coordinate system can be used as input and output.
%
%   PLH=MAP2PLH(MAP,CRSIN,CRSOUT,METHOD) allows the user to specify
%   the transformation tool. Current options are 'proj4' and 
%   'rdnaptrans' (latter for transformation to and from the Dutch
%   RDNAP system only). Default is 'proj4';
%
%   Examples:                                   
%       plh = crstrans(xyh,'RD','WGS84')
%       plh = crstrans(xyh,'EPGS28992','WGS84','rdnaptrans')
%       xyh = crstrans(plh,'WGS84','EPGS28992','rdnaptrans')
%     
%   See also CS2CS, RDNAP2ETRS, ETRS2RDNAP.
%
%   (c) Freek van Leijen, Delft University of Technology, 2022.

%   Created:    20 Jan 2022 by Freek van Leijen
%   Modified:   

% Input argument checking

if nargin < 3
  error('Must be called with at least three arguments.');
end

if nargin < 4
  method = 'proj4';
end

if size(crdIn,2) ==2
  crdIn = [crdIn zeros(size(crdIn,1),1)];
elseif size(crdIn,2) ~= 3
  error('CRDIN must have two or three columns or elements.')
end

crsIn = upper(crsIn);
crsOut = upper(crsOut);


%% Transformation

switch lower(method)
  case 'proj4'
    [crd1,crd2] = cs2cs(crsIn,crsOut,crdIn(:,1),crdIn(:,2));
    crdOut = [crd1 crd2 zeros(size(crd1))];
  case 'rdnaptrans'
    if strmatch(crsIn,{'RD','RDNAP','EPSG28992'}) & ...
       strmatch(crsOut,{'WGS84','ETRS','ITRS','EPSG4326'})
      crdOut = rdnap2etrs(crdIn,'PLH');
      crdOut(:,1:2) = crdOut(:,1:2)*180/pi;
    elseif strmatch(crsIn,{'WGS84','ETRS','ITRS','EPSG4326'}) & ...
       strmatch(crsOut,{'RD','RDNAP','EPSG28992'})
      crdOut = etrs2rdnap([crdIn(:,1:2)*pi/180 crdIn(:,3)],'PLH');
    else
      error('Not possible to perform this transformation with RDNAPTRANS.');
    end
  otherwise
    error('You specified an unsupported coordinate transformation tool.');
end

end
