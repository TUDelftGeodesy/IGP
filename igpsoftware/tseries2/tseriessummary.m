function ht=tseriessummary(stations,id,ptitle)
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
%   HT=TSERIESSUMMARY(STATIONS,ID,PTITLE) in addition to printing a summary
%   of the estimated time series parameter, also makes two bar plots
%   of the relevant parameters, with plot title PTITLE. The output HT
%   is an array with the plot handles.
%
%   Examples:
%      tseriessummary(tseries);
%      tseriessummary({'delf' 'kosg' 'wsrt' 'ters' 'eijs'});
%      tseriessummary(stations,'_fit2');
%      tseriessummary(stations,'_fit2','Plot title');
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016, 2021.

% Check arguments

if nargin < 1, error('too few arguments (specify at least stationames');end
if nargin > 3, error('too many arguments');end
if ~iscell(stations) && ~isstruct(stations)
    error('first argument must be a cell array with station names or a structure array of timeseries.') 
end
if nargin < 2
    id='_fit';
end
if nargin > 2
  doplot=true;
else
  doplot=false;
end

% number of stations

nstations=numel(stations);          

% create table legend

xlabel={ 'Vel1 [mm/y]' 'Vel2 [mm/y]' 'AtmLd [mm/kPa]' ,'Temp. [mm/daK]' } ;

xlabel1='             Vel1   Vel2     AtmLd  TempI  ';
xlabel2='             mm/y   mm/y    mm/kPa mm/daK  ';

maxharmonic=0;
for i=1:nstations
  if iscell(stations)
    station=stations{i};
    load([ station id '.mat']);
  elseif isstruct(stations)
    tseries=stations(i);
    station=tseries.station;
  end
  nharmonic=length(tseries.harmonic);
  for k=1:nharmonic
     xlabel{4+k}=sprintf('%3dd [mm]',round(tseries.harmonic(k)*365));
     xlabel3{k}=sprintf('   %3dd',round(tseries.harmonic(k)*365));
   end
  maxharmonic=max(maxharmonic,nharmonic);
end

fmt='';
for k=1:maxharmonic
   xlabel1=[ xlabel1 xlabel3{k} ];
   xlabel2=[ xlabel2 '     mm' ];
   fmt=[ fmt ' %6.2f'];
end

xlabel=[ xlabel {'StdF [mm]' 'StdE [mm]' 'StdR [mm]' 'OMT [-]' }];
xlabel1= [ xlabel1 '     StdF   StdE   StdR      OMT' ];
xlabel2= [ xlabel2 '       mm     mm     mm         ' ];

fid=1;
fprintf(fid,'%s\n',xlabel1);
fprintf(fid,'%s\n\n',xlabel2);

% Loop over stations and print the parameteres

nmax=8+nharmonic;

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
  
  stationlist{i}=station;
  
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

  nharmonic=length(tseries.harmonic);

  sumgen(i,1:2) = [ tseries.tfirst tseries.tlast ];
  sumlat(i,1:2) = [ tseries.vitrf(1)*1000 tseries.vel(1)*1000 ];
  sumlon(i,1:2) = [ tseries.vitrf(2)*1000 tseries.vel(2)*1000 ];
  sumrad(i,1:2) = [ tseries.vitrf(3)*1000 tseries.vel(3)*1000 ];
  nmeteo=parindex(2)-parindex(1);
  for j=1:nmeteo
     sumlat(i,2+j) = tseries.xlat(npar+j)*1000;
     sumlon(i,2+j) = tseries.xlon(npar+j)*1000;
     sumrad(i,2+j) = tseries.xrad(npar+j)*1000;
  end
  for j=1:nharmonic
     jj=4+j;
     sumlat(i,jj) = tseries.amplitude(1,j)*1000;
     sumlon(i,jj) = tseries.amplitude(2,j)*1000;
     sumrad(i,jj) = tseries.amplitude(3,j)*1000;
  end
  sumlat(i,maxharmonic+5:maxharmonic+8)= [ tseries.meanstd(1)*1000 tseries.empstd(1)*1000 tseries.rms(1)*1000 tseries.omt(1) ];  
  sumlon(i,maxharmonic+5:maxharmonic+8)= [ tseries.meanstd(2)*1000 tseries.empstd(2)*1000 tseries.rms(2)*1000 tseries.omt(2) ];  
  sumrad(i,maxharmonic+5:maxharmonic+8)= [ tseries.meanstd(3)*1000 tseries.empstd(3)*1000 tseries.rms(3)*1000 tseries.omt(3) ];  

  fprintf(fid,[ '%s  Lat  %6.2f %6.2f    %6.2f %6.2f  ' fmt '   %6.2f %6.2f %6.2f   %6.2f\n'],station,sumlat(i,1:nmax));
  fprintf(fid,[ '%s  Lon  %6.2f %6.2f    %6.2f %6.2f  ' fmt '   %6.2f %6.2f %6.2f   %6.2f\n'],station,sumlon(i,1:nmax));
  fprintf(fid,[ '%s  Rad  %6.2f %6.2f    %6.2f %6.2f  ' fmt '   %6.2f %6.2f %6.2f   %6.2f\n'],station,sumrad(i,1:nmax));
  fprintf(fid,'\n');

end

% Optional bar charts  

if doplot

    pid=id;
    pid2=id;
    ptitle2=ptitle;
    if ~isempty(id)
      pid=[ '_' id];
      pid2=[' (it=' id ') - '];
    else
      ptitle2=[ ptitle2 ' - ' ];
    end

    nsel= [ 3:8 nmax-1 nmax ];
    for k=1:8
      ksel=nsel(k);
      %  ax=subplot(1,nsel,ksel);
      h(k)=figure;
      barh([sumlat(:,ksel) sumlon(:,ksel) sumrad(:,ksel)]);
      ylim([0 nstations+1])
      set(gca,'YTick',1:nstations,'YTickLabel',stationlist,'YDir','reverse')
      ax(k)=gca;
      title([ ptitle2 pid2 xlabel{ksel} ]);
    end


    ht(1)=fig2subplot([ h(1) h(3) ; h(2) h(4)]);
    legend({'North' 'East' 'Up'})
    set(ht(1),'Name',[ ptitle '_bar_parms' pid ],'NumberTitle','off');

    ht(2)=fig2subplot([ h(5) h(6) ; h(7) h(8)]);
    legend({'North' 'East' 'Up'})
    set(ht(2),'Name',[ ptitle '_bar_stats' pid ],'NumberTitle','off');

    close(h);

end


end

