function tseriessummary(stations,id)
%TSERIESSUMMARY Print summary of estimated parameters and standard deviations.
%   TSERIESSUMMARY(TSERIES) print a summary of the estimated time
%   series parameters for the stations in the structure array TSERIES.
%
%   TSERIESSUMMARY(STATIONS,ID) print a summary of the estimated time
%   series parameters for the stations in the cell array STATIONS. ID is 
%   the optional identifier of the fit (default '_fit'). The time series
%   data must be stored in an mat-file STATION_ID.mat (_ is part of the
%   ID).
%
%   Examples:
%      tseriessummary(tseries);
%      tseriessummary({'delf' 'kosg' 'wsrt' 'ters' 'eijs'});
%      tseriessummary(stations,'_fit2');
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016.

% Check arguments

if nargin < 1, error('too few arguments (specify at least stationames');end
if nargin > 2, error('too many arguments');end
if ~iscell(stations) && ~isstruct(stations)
    error('first argument must be a cell array with station names or a structure array of timeseries.') 
end
if nargin < 2
    id='_fit';
end

% Loop over stations and print the parameteres

nstations=numel(stations);          

nmax=40;

sumgen=zeros(nstations,2);
sumlat=zeros(nstations,nmax);
sumlon=zeros(nstations,nmax);
sumrad=zeros(nstations,nmax);

for i=1:nstations

  if iscell(stations)
    station=stations{i};
    load([ station id '.mat']);
  elseif isstruct(stations)
    tseries=stations(i);
    station=tseries.station;
  end
  parlabel=tseries.parlabel;
  parindex=tseries.parindex;

  % correct some mistakes in the input structure (now corrected)

  if ~isfield(tseries,'virtf'), tseries.vitrf=[ 0 0 0];, end
  if ~isfield(tseries,'empstd') || ~isfield(tseries,'meanstd')
    tseries.empstd=std(diff(tseries.neuresidual))./sqrt(2);
    tseries.meanstd=mean(tseries.sneu);
  end
  
  % special treatment for splines
  
  npar=parindex(1);
  if isfield(tseries,'knots') && ~isempty(tseries.knots)
    spline=true;
  else
    spline=false;
  end
  if spline
    sp=spmak(tseries.knots,[tseries.xlat(1:npar)'; tseries.xlon(1:npar)'; tseries.xrad(1:npar)']);
    spvel=fnder(sp);  
    tseries.vel=fnval(spvel,0)';
  end

  % create table legend
  
  xlabel={ 'Vel1 [mm/y]' 'Vel2 [mm/y]' 'AtmLd [mm/kPa]' ,'Temp. [mm/daK]' } ;
  xlabel1='             Vel1   Vel2     AtmLd  TempI  ';
  xlabel2='             mm/y   mm/y    mm/kPa mm/daK  ';
  nharmonic=length(tseries.harmonic);
  fmt='';
  for k=1:nharmonic
     xlabel{4+k}=sprintf('%3dd [mm]',round(tseries.harmonic(k)*365));
     xlabel1=[ xlabel1 sprintf('   %3dd',round(tseries.harmonic(k)*365))];
     xlabel2=[ xlabel2 '     mm' ];
     fmt=[ fmt ' %6.2f'];
  end
  xlabel=[ xlabel {'StdF [mm]' 'StdE [mm]' 'StdR [mm]' 'OMT [-]' }];
  xlabel1= [ xlabel1 '     StdF   StdE   StdR      OMT' ];
  xlabel2= [ xlabel2 '       mm     mm     mm         ' ];

  npar=parindex(1);
  n=nharmonic+8;
  sumgen(i,1:2) = [ tseries.tfirst tseries.tlast ];
  sumlat(i,1:n) = [ tseries.vitrf(1)*1000 tseries.vel(1)*1000 tseries.xlat(npar+1)*1000 tseries.xlat(npar+2)*1000 tseries.amplitude(1,:)*1000 tseries.meanstd(1)*1000 tseries.empstd(1)*1000 tseries.rms(1)*1000 tseries.omt(1) ];  
  sumlon(i,1:n) = [ tseries.vitrf(2)*1000 tseries.vel(2)*1000 tseries.xlon(npar+1)*1000 tseries.xlon(npar+2)*1000 tseries.amplitude(2,:)*1000 tseries.meanstd(2)*1000 tseries.empstd(2)*1000 tseries.rms(2)*1000 tseries.omt(2) ];  
  sumrad(i,1:n) = [ tseries.vitrf(3)*1000 tseries.vel(3)*1000 tseries.xrad(npar+1)*1000 tseries.xrad(npar+2)*1000 tseries.amplitude(3,:)*1000 tseries.meanstd(3)*1000 tseries.empstd(3)*1000 tseries.rms(3)*1000 tseries.omt(3) ];  

  fid=1;
  if i==1 
    fprintf(fid,'%s\n',xlabel1);
    fprintf(fid,'%s\n\n',xlabel2);
  end
  fprintf(fid,[ '%s  Lat  %6.2f %6.2f    %6.2f %6.2f  ' fmt '   %6.2f %6.2f %6.2f   %6.2f\n'],station,sumlat(i,1:n));
  fprintf(fid,[ '%s  Lon  %6.2f %6.2f    %6.2f %6.2f  ' fmt '   %6.2f %6.2f %6.2f   %6.2f\n'],station,sumlon(i,1:n));
  fprintf(fid,[ '%s  Rad  %6.2f %6.2f    %6.2f %6.2f  ' fmt '   %6.2f %6.2f %6.2f   %6.2f\n'],station,sumrad(i,1:n));
  fprintf(fid,'\n');
  %sumlat(i,1:n)
  %sumlon(i,1:n)
  %sumrad(i,1:n)

end

end

