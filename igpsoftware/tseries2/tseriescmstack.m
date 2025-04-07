function tseries=tseriescmstack(rstack,cm,env)
%tseriesstack    Compute timeseries structure from the stack data.
%   TS=TSERIESCMSTACK(RSTACK,CM,ENV) builds a timeseries structure TS from
%   the residual stack in the RSTACK structure and common mode in the
%   stucture CM, for the epochs in the RSTACK structure. ENV is an optionalt
%   two dimensional array with the temperature and pressure. Output is an 
%   time series structuere TS.
%
%   TS=TSERIESCMSTACK(RSTACK,CM,METSTRUCT) does the same, but retrieved the
%   meteodata from a metel structure METSTRUCT or meteo file with name
%   METSTRUCT.
%
%   TS=TSERIESCMSTACK(RSTACK,CM) evaluates the common mode only for the
%   harmonic components.
%
%   Examples:
%      ts=tseriescmstack(rstack,cm)
%      ts=tseriescmstack(rstack,cm,metstruct);
%      ts=tseriescmstack(rstack,cm,'meteo_Eelde.mat');
%
%   See also TSERIESRESIDUALSTACK, TSERIESCMFIT and TSERIESCMEVAL.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2019.

% Check arguments

if nargin < 3; env=[]; end

% Check if residual stack needs loading

if ischar(rstack)
  load(rstack);
end

titlestr=rstack.titlestr;
stations=rstack.stations;

year=rstack.ryear;
rmean=rstack.rmean;
rstd=rstack.rstd;
rcount=rstack.rcount;

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
   for k=1:numel(cm.parlabel)
      cm.parunit{k}='mm';
   end
   if isempty(iAtmLd), iharmonic=iharmonic-1; else cm.parunit{iAtmLd}='mm/kPa'; end
   if isempty(iTempI), iharmonic=iharmonic-1; else cm.parunit{iAtmLd}='mm/daK'; end
   cm.t0=round(mean(year));
   cm.parindex(1)=0;
   cm.parindex(2)=iharmonic-1;
   cm.parindex(3)=numel(cm.parlabel);
   cm.parindex(4)=cm.parindex(3);
end

% Save results in a structure

tseries.station='CM';
tseries.id='';
tseries.plh0=[nan nan nan];
tseries.vneu0=[nan nan nan];
tseries.events=[];

tseries.tfirst=min(year);
tseries.tlast=max(year);
tseries.t0=cm.t0;

tseries.year=year;
tseries.neu=rmean;
tseries.sneu=rstd;

tseries.model='zero-velocity';
tseries.method='polynomial';
tseries.pdegree=-1;
tseries.env=env;

tseries.env0=cm.env0;
tseries.harmonic=cm.harmonic;
tseries.tjump=[];

tseries.parindex=cm.parindex;
tseries.parlabel=cm.parlabel;
tseries.parunit=cm.parunit;

tseries.xlat=cm.xlat;
tseries.xlon=cm.xlon;
tseries.xrad=cm.xrad;

tseries.qxlat=cm.qxlat;
tseries.qxlon=cm.qxlon;
tseries.qxrad=cm.qxrad;

tseries.meanstd=[nan nan nan ];
tseries.empstd=[nan nan nan ];
tseries.rms=[nan nan nan ];
tseries.omt=[nan nan nan ];

tseries.knots=[];
tseries.pspline=[];
tseries.sp=[];
tseries.spvel=[];
tseries.vel=[ 0 0 0 ];
tseries.svel=[ 0 0 0 ];

tseries.lattext='';
tseries.lontext='';
tseries.radtext='';

amplitude=zeros(3,length(cm.harmonic));
npar=cm.parindex(2);
for i=1:length(cm.harmonic)
   amplitude(1,i)=sqrt(cm.xlat(npar+2*(i-1)+1)^2+cm.xlat(npar+2*(i-1)+2)^2);
   amplitude(2,i)=sqrt(cm.xlon(npar+2*(i-1)+1)^2+cm.xlon(npar+2*(i-1)+2)^2);
   amplitude(3,i)=sqrt(cm.xrad(npar+2*(i-1)+1)^2+cm.xrad(npar+2*(i-1)+2)^2);
end
tseries.amplitude=amplitude;

% Get the temperature and pressure data

if isempty(env)
  P=ones(size(year))*cm.env0(1);   
  T=ones(size(year))*cm.env0(2);
elseif isstruct(env) || ischar(env);
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
  P=env{1};
  T=env{2};
end

% Time series residuals (residual stack)

tseries.neuresidual=rmean;    

% Zero trend 

tseries.neutrend=zeros(size(rmean));

% Temperature and atmospheric loading time series

m=size(year,1);
if ~isempty(P) && ~isempty(iAtmLd)
  dP=0.1*(P-cm.env0(1));
  tseries.neuatmld=  [ dP * cm.xlat(iAtmLd)  ...
                       dP * cm.xlon(iAtmLd)  ...
                       dP * cm.xrad(iAtmLd) ];
else
  tseries.neuatmld=zeros(m,3);
end
if ~isempty(T) && ~isempty(iTempI)
  dT=0.1*(T-cm.env0(2));
  tseries.neutempi=  [ dT * cm.xlat(iTempI)  ...
                       dT * cm.xlon(iTempI)  ...
                       dT * cm.xrad(iTempI) ];
else
  tseries.neutempi=zeros(m,3);
end

% Harmonic time series

t0=cm.t0;
dt=year-t0;
nenv=cm.parindex(2);
tseries.neuharmonic=zeros(m,3);
for i=1:length(cm.harmonic)
   f=1/cm.harmonic(i);
   tseries.neuharmonic=tseries.neuharmonic + ...
       [  sin(2*pi*f*dt)*cm.xlat(nenv+1)+cos(2*pi*f*dt)*cm.xlat(nenv+2) ...
          sin(2*pi*f*dt)*cm.xlon(nenv+1)+cos(2*pi*f*dt)*cm.xlon(nenv+2) ...
          sin(2*pi*f*dt)*cm.xrad(nenv+1)+cos(2*pi*f*dt)*cm.xrad(nenv+2) ];
   nenv=nenv+2;
end

% No jumps in the common mode

tseries.neujumps=zeros(size(rmean));

% Fit

tseries.neufit=tseries.neutrend+tseries.neuatmld+tseries.neutempi+tseries.neuharmonic;


end





