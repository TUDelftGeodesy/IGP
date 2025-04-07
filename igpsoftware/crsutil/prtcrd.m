function prtcrd(station,crd,crdformat,varargin)
%PRTCRD  Print a table with coordinates and optional co-variance information.
%   PRTCRD(STATION,CRD) prints a table with the station names STATION and 
%   coordinates given by the M-by-3 matrix CRD. 
%
%   PRTCRD(STATION,CRD,STDEV) prints three extra columns with the standard
%   deviations from the M-by-3 matrix STDEV.
%
%   PRTCRD(STATION,CRD,CRDFMT) uses the format CRDFMT for printing.
%   The coordinates in CRD are assumed to match CRDFMT. Supported values 
%   for CRDFMT are
%
%      'ecef' 'xyz'        Cartesian XYZ coordinates   (default)
%      'local' 'neu'       Coordinates in local NEU system
%      'map'               Coordinates in map projection system
%      'geodetic' 'plh'    Geographic coordinates printed in decimal degrees
%      'sexagesimal' 'dms' Geographic coordiantes printed in degrees, minutes, 
%                          seconds
%
%   For the options 'dms' and 'plh' CRD must contain the latitude and 
%   longitude in radians.
%
%   PRTCRD(STATION,CRD,SDCOV,CRDFMT,COVFMT) uses the format COVFMT for 
%   of printing the standard deviations and co-variances in the M-by-3
%   or M-by-6 matrix SDCOV. Valid formats for COVFMT are 'std' to print
%   only standard deviations, or any format supported by COVREFORMAT
%   function (except 'qmat').
%
%   PRTCRD(...,TITLE) with TITLE as 4th or 6th argument print an title
%   above the table.
%
%   See also COVREFORMAT and STR2PLH.
%
%   (c) Hans van der Marel, Delft University of Technology, 2014

%   Created:    12 March 2014 by Hans van der Marel
%   Modified:   

% Defaults

sdcov=[];
covformat=[];
title=[];

% Checking of input arguments and options

if nargin < 2
   error('insufficient number of arguments')
elseif nargin == 2 
   crdformat='xyz';
elseif nargin == 3 || nargin == 4
   if ~ischar(crdformat)
       sdcov=crdformat;
       crdformat='xyz';
       covformat='std';
   end      
   if nargin == 4
     title=varargin{1};
   end
elseif nargin == 5 || nargin == 6
   sdcov=crdformat;
   crdformat=varargin{1};
   covformat=varargin{2};
   if nargin == 6
       title=varargin{3};
   end
else
  error('incorrect number of arguments')
end

if size(crd,2) ~= 3
  error('crd array must have three columns with coordinates')
end
if ~isempty(sdcov) 
  if size(sdcov,2) == 3
     covformat='std';
  elseif size(sdcov,2) ~= 6
     error('sdcov array must have three or six columns')    
  end
end

% Convert station names into character array

if iscell(station)            
  station=char(station);      
elseif ~ischar(station)
  station=num2str(station,'%d');
end

% Number of stations

nsta=size(crd,1);            
if size(station,1) ~= nsta
  error('size of input arrays does not match (station and crd)')  
end
if ~isempty(sdcov) && size(sdcov,1) ~= nsta
  error('size of input arrays does not match (crd and sdcov)')  
end

% Format coordinates

switch lower(crdformat)       
    case {'xyz','ecef'}
       crdheader='            X[m]            Y[m]            Z[m]';
       covheader='    sx[m]    sy[m]    sz[m]   sxy[m]   sxz[m]   syz[m]';
       crdstr=reshape(sprintf('%16.4f',crd'),[48,nsta])';
    case {'neu','local'}
       crdheader='        N[m]        E[m]        U[m]';
       covheader='    sn[m]    se[m]    su[m]   sne[m]   snu[m]   seu[m]';
       crdstr=reshape(sprintf('%12.4f',crd'),[36,nsta])';
    case {'map'}
       crdheader='          x[m]          y[m]        h[m]';
       covheader='    sx[m]    sy[m]    sh[m]   sxy[m]   sxh[m]   syh[m]';
       crdstr=reshape(sprintf('%14.4f%14.4f%12.4f',crd'),[40,nsta])';
    case {'plh','geodetic'}
       crdheader='        Lat[deg]        Lon[deg]       H[m]';
       covheader='    sn[m]    se[m]    su[m]   sne[m]   snu[m]   seu[m]';
       crdstr=reshape(sprintf('%16.8f%16.8f%11.4f',[ crd(:,1:2)*180/pi crd(:,3) ]'),[43,nsta])';
    case {'dms','sexagesimal'}
       crdheader='        Lat[dms]        Lon[dms]       H[m]';
       covheader='    sn[m]    se[m]    su[m]   sne[m]   snu[m]   seu[m]';
       crdstr=plh2str(crd);
    otherwise
       error('unknown coordinate format')
end

% Format covariance information

if ~isempty(sdcov)            
    switch covformat
        case 'std'
           covheader=covheader(1:27);
           covstr=reshape(sprintf('%9.4f',sdcov(:,1:3)'),[27,nsta])';
        otherwise
           sdcov=covreformat(sdcov,covformat,'scov');
           covstr=reshape(sprintf('%9.4f',sdcov'),[54,nsta])';
    end
else
    covheader=[];
    covstr=[];
end

% Print

if ~isempty(title)
  fprintf('%s\n',title) 
end
disp([ repmat(' ',[1,size(station,2)]) crdheader '   ' covheader ])
disp([ station  crdstr repmat(' ',[nsta,3]) covstr ])
fprintf('\n')

return