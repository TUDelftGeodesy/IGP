function y=tseriescmeval(cm,year,env)
%TSERIESCMEVAL  Evaluate common mode fit.
%   Y=TSERIESCMEVAL(CM,YEAR,ENV) evaluates the common mode in the
%   stucture CM for the decimal year array YEAR. ENV is a two dimensional
%   array with the temperature and pressure. Output is an array Y with
%   three columns for the common mode in the latitude, longitude and 
%   vertical component. The input CM is a structure with the estimated
%   common mode.
%
%   Y=TSERIESCMEVAL(CM,YEAR,METSTRUCT) does the same, but retrieved the
%   meteodata from a metel structure METSTRUCT or meteo file with name
%   METSTRUCT.
%
%   Y=TSERIESCMEVAL(CM,YEAR) evaluates the common mode only for the
%   harmonic components.
%
%   Y=TSERIESCMEVAL(TSERIES,...) evaluates the loading and harmonic effects
%   for the timeseries TSERIES.
%
%   TSERIESCMEVAL(...) plots the common mode.
%
%   Examples:
%      cm=tseriescmfit(stations,,'_fit');
%      y=tseriescmeval(cm,year,env);
%      tseriescmeval(cm,year,env);
%      tseriescmeval(tseries,year,'meteo_Eelde.mat');
%
%   See also TSERIESCMFIT.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2019.

% Check arguments

if nargin < 3; env=[]; end

alat=[];
alon=[];
arad=[];

year=year(:);

% Check if we have a cm or tseries structure

iAtmLd=find(strcmpi(cm.parlabel,'AtmLd'));
iTempI=find(strcmpi(cm.parlabel,'TempI'));
if isfield(cm,'parindex')
   % tseries structure
   iharmonic=cm.parindex(2)+1; 
   station=cm.station;
   if ~isfield(cm,'env0')   % backward compability 
     cm.env0=[ mean(cm.env{1}) mean(cm.env{2})];
   end
else
   % cm stucture 
   iharmonic=3;  
   if isempty(iAtmLd), iharmonic=iharmonic-1;, end
   if isempty(iTempI), iharmonic=iharmonic-1;, end
   station='Common Mode';
   cm.t0=round(mean(year));
end

% Get the temperature and pressure data

if isempty(env)
  subtitle='Harmonics';
  P=ones(size(year))*cm.env0(1);   
  T=ones(size(year))*cm.env0(2);
elseif isstruct(env) || ischar(env);
  subtitle='Harmonics+Env';
  % meteo structure
  if ischar(env)
    metstruct=load(env);
  else
    metstruct=env;
  end
  YMD=dyear2ymd(year);
  d2=datestr(datenum(YMD),'yyyymmdd');
  [~,k2]=ismember(d2,metstruct.YYYYMMDD,'rows');
  P=metstruct.Pday(k2);
  T=metstruct.Tday(k2);
else
  subtitle='Harmonics+Env';
  P=env{1};
  T=env{2};
end

% Set up the temperature and atmospheric loading model

xlat=[];xlon=[];xrad=[];
if ~isempty(P) && ~isempty(iAtmLd)
  dP=0.1*(P-cm.env0(1));
  alat = [alat dP ];
  alon = [alon dP ];
  arad = [arad dP ];
  xlat=[xlat; cm.xlat(iAtmLd)];
  xlon=[xlon; cm.xlon(iAtmLd)];
  xrad=[xrad; cm.xrad(iAtmLd)];
end
if ~isempty(T) && ~isempty(iTempI)
  dT=0.1*(T-cm.env0(2));
  alat = [alat dT ];
  alon = [alon dT ];
  arad = [arad dT ];
  xlat=[xlat; cm.xlat(iTempI)];
  xlon=[xlon; cm.xlon(iTempI)];
  xrad=[xrad; cm.xrad(iTempI)];
end
xlat=[xlat; cm.xlat(iharmonic:iharmonic+length(cm.harmonic)*2-1)];
xlon=[xlon; cm.xlon(iharmonic:iharmonic+length(cm.harmonic)*2-1)];
xrad=[xrad; cm.xrad(iharmonic:iharmonic+length(cm.harmonic)*2-1)];

% Set up the harmonics model

%t0=round(mean(year));
t0=cm.t0;
dt=year-t0;
for i=1:length(cm.harmonic)
   f=1/cm.harmonic(i);
   alat=[alat sin(2*pi*f*dt) cos(2*pi*f*dt) ];
   alon=[alon sin(2*pi*f*dt) cos(2*pi*f*dt) ];
   arad=[arad sin(2*pi*f*dt) cos(2*pi*f*dt) ];
end

% Evaluate

y=[ alat*xlat alon*xlon arad*xrad ];

% Plot if no output arguments given

if nargout < 1
  figure
  plot(year,y*1000)
  ylabel('[mm]')
  title(['Periodic effects: ' station ' (' subtitle ')'])
  legend('North','East','Up')
end

end
