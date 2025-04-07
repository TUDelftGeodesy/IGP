function strvalue=proj4defs(parameter)
%PROJ4DEFS   Return default values for PROJ.4 projection support
%   STRVALUE=PROJ4DEFS(PARAMETER) returns the value STRVALUE for PROJ.4
%   projection library settings. PARAMETER is a text string with the
%   requested parameter:
%
%      PROJID       returns the PROJ.4 +args for the projection PROJID (e.g.
%                   'RD', 'WGS84', ...
%      EPSG######   returns the PROJ.4 +args for the projectcion with
%                   EPSG code ######
%      'PROJ4PATH'  returns path to PROJ.4 library and PROJ and CS2CS programs
%
%   STRVALUE=PROJ4DEFS() returns a cell array with projection names
%   for PROJ.4. 
%
%   PROJ4DEFS() prints projection names for PROJ.4. 
%
%   Example:
%
%       proj4path=proj4defs('PROJ4PATH')
%       proj4pars=proj4defs('RD')
%       proj4pars=proj4defs('EPGS28992')
%       proj4defs()
%
%   See also proj4fwd, proj4inv and cs2cs.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015.
%
%   Modified:  21 September 2023 by Freek van Leijen
%               - changed EPSG4326 proj definition from longlat to latlong

% Definitions

proj.WGS84='EPSG4326';
proj.RD='EPSG28992';
proj.DHDN2='EPSG31466';
proj.LB72='EPSG31700';   % 'EPSG31300' 'EPSG103300'

% WGS84
%proj.EPSG4326='+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs';
proj.EPSG4326='+proj=latlong +ellps=WGS84 +datum=WGS84 +no_defs';

% Amersfoort / RD New
proj.EPSG28992='+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.417,50.3319,465.552,-0.398957,0.343988,-1.8774,4.0725 +units=m +no_defs';       % Amersfoort / RD New    (From cs2cs webtool, older epgs dist)
%roj.EPSG28992='+proj=sterea +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.9999079 +x_0=155000 +y_0=463000 +ellps=bessel +towgs84=565.4171,50.3319,465.5524,-0.398957388243134,0.343987817378283,-1.87740163998045,4.0725 +units=m +no_defs';    % Amersfoort / RD New    (From epsg 4.9.1 dist)
%roj.EPSG28992='+proj=stere +lat_0=52.15616055555555 +lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel +units=m +no_defs  no_defs';    % # Amersfoort / RD New    (From esri 4.9.1 dist) 

% DHDN / 3-degree Gauss-Kruger zone 2
proj.EPSG31466='+proj=tmerc +lat_0=0 +lon_0=6 +k=1 +x_0=2500000 +y_0=0 +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7 +units=m +no_defs';  % DHDN / 3-degree Gauss-Kruger zone 2 (From cs2cs webtool)
%roj.EPSG31466='+proj=tmerc +lat_0=0 +lon_0=6 +k=1 +x_0=2500000 +y_0=0 +datum=potsdam +units=m +no_defs';   % DHDN / 3-degree Gauss-Kruger zone 2 (From esri 4.9.1 dist) +datum=potsdam --> +ellps=bessel +towgs84=598.1,73.7,418.2,0.202,0.045,-2.455,6.7

% Belge Lambert 72
proj.EPSG31300='+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=90 +lon_0=4.356939722222222 +x_0=150000.01256 +y_0=5400088.4378 +ellps=intl +towgs84=-106.868628,52.297783,-103.723893,0.336570,-0.456955,1.842183,-1.2747 +units=m +no_defs';  % Belge 1972 / Belge Lambert 72  (From cs2cs webtool)
proj.EPSG31700='+proj=lcc +lat_1=51.16666723333333 +lat_2=49.8333339 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.013 +y_0=5400088.438  +ellps=intl +towgs84=-106.868628,52.297783,-103.723893,0.336570,-0.456955,1.842183,-1.2747 +units=m +no_defs';  % Belge 1972 / Belge Lambert 72  (From cs2cs webtool)
%roj.EPSG31300='+proj=lcc +lat_1=49.83333333333334 +lat_2=51.16666666666666 +lat_0=90 +lon_0=4.356939722222222 +x_0=150000.01256 +y_0=5400088.4378 +ellps=intl +towgs84=-106.8686,52.2978,-103.7239,0.3366,-0.457,1.8422,-1.2747 +units=m +no_defs'; % Belge 1972 / Belge Lambert 72  (From epsg file)
%roj.EPSG31700='+proj=lcc +lat_1=51.16666723333333 +lat_2=49.8333339 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.013 +y_0=5400088.438 +ellps=intl +towgs84=-106.8686,52.2978,-103.7239,0.3366,-0.457,1.8422,-1.2747 +units=m +no_defs'; % Belge 1972 / Belge Lambert 72  (From epsg file)
proj.EPSG103300='+proj=lcc +lat_1=49.8333339 +lat_2=51.16666733333333 +lat_0=90 +lon_0=4.367486666666666 +x_0=150000.01256 +y_0=5400088.4378 +ellps=intl +units=m +no_defs'; 


% check input arguments

if nargin > 1
  error('function needs one string parameter or nothing')
end

% print list of projections

if nargin < 1
  projnames=fieldnames(proj);
  if nargout > 0
    strvalue=projnames;
  else
    for k=1:length(projnames)
      fprintf('%s\n',projnames{k})
    end
  end
  return
end

% process input argument

ucparameter=upper(parameter);
switch ucparameter
  case 'PROJ4PATH'
    thisdir=fileparts(mfilename('fullpath'));
    if exist(fullfile(thisdir,'proj4path.mat'),'file')
       load(fullfile(thisdir,'proj4path'));
    else
       % try to find PROJ.4 (e.g. QGIS install) - Windows only
       searchloc= { fullfile(getenv('ProgramFiles'),'QGIS*') ; ...
                    fullfile(getenv('ProgramW6432'),'QGIS*')};
       optionpath={};ksel=0;
       for l=1:length(searchloc)
         d=dir(searchloc{l});
         for k=1:length(d)
            ksel=ksel+1;
            optionpath{ksel}=fullfile(fileparts(searchloc{l}),d(k).name,'bin');
         end
       end
       if ksel > 1
         fprintf('Possible locations of PROJ.4:\n')
         for k=1:ksel
            fprintf('%d: %s\n%',k, optionpath{k});
         end
         ksel=input('Your choice:');
       end
       if ksel >= 1
          PROJ4PATH=optionpath{ksel};
       else 
         PROJ4PATH=input('Please enter path to directory with PROJ.4 executable:','s');
       end
       if exist(PROJ4PATH,'dir')
          fprintf('Set PROJ.4 path to %s\n',PROJ4PATH)
          save(fullfile(thisdir,'proj4path'),'PROJ4PATH')
       end
    end
    if ~exist(PROJ4PATH,'dir')
       error('Path to PROJ.4 executables could not be determined')
    end
    %PROJ4PATH='c:\Program Files\QGIS Chugiak\bin\';
    strvalue=PROJ4PATH;
  case 'PROJEXE'
    PROJ4PATH=proj4defs('PROJ4PATH');
    if ispc
      strvalue=fullfile(PROJ4PATH,'proj.exe');
    else
      strvalue=fullfile(PROJ4PATH,'proj');
    end
  case 'CS2CSEXE'
    PROJ4PATH=proj4defs('PROJ4PATH');
    if ispc
      strvalue=fullfile(PROJ4PATH,'cs2cs.exe');
    else
      strvalue=fullfile(PROJ4PATH,'cs2cs');
    end
  otherwise
    strvalue=parameter;
    while ~strncmp(parameter,'+',1)
      if isfield(proj,ucparameter)
        strvalue=proj.(ucparameter);
        parameter=strvalue;
        ucparameter=upper(parameter);
      else
        error('illegal parameter')
      end
    end
end

end


