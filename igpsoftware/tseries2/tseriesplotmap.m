function hf=tseriesplotmap(stations,item,id,project)
%TSERIESPLOTMAP  Plot map with representation of the estimated parameters.
%   TSERIESPLOTMAP(TSERIES,ITEM,ID) plot map with the timeseries components
%   in the cell array ITEM for the stations in the structure array TSERIES.
%   ID is an optional identifier.
%
%   TSERIESPLOTMAP(STATIONS,ITEM,ID) plot map with the timeseries components
%   in the cell array ITEM for the stations in the cell array STATIONS. 
%   ID is the optional identifier of the fit (default '_fit'). The timeseries
%   data must be stored in an mat-file STATION_ID.mat (_ is part of the
%   ID).
%
%   TSERIESPLOTMAP(...,PROJECT) provides additional project
%   name PROJECT.
%
%   Examples:
%      tseriesplotmap(tseries,'raw');
%      tseriesplotmap(stations,{'harmonic' 'tempi'},'_fit');
%      tseriesplotmap(stations,{'harmonic' 'tempi'},'_fit','Ameland');
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016, 2021.

% Check arguments

if nargin < 1, error('too few arguments (specify at least stationames and item');end
if ~iscell(stations) && ~isstruct(stations)
    error('first argument must be a cell array with station names or a structure array of timeseries.') 
end
if nargin < 2
    item='maponly';
end
if nargin < 3
    id='_fit';
end
if nargin < 4
   project='';
end

saveplots=false;

titlestr=[project ' Map'];

nstations=length(stations);          

if iscell(stations)
   stationnames=stations;
   plh0=zeros(nstations,3);
   for i=1:nstations
     %station=stations{i};
     %load([ station '.mat']);
     station=stations{i};
     load([ station id '.mat']);
     plh0(i,:)=tseries.plh0;
   end
elseif isstruct(stations)
   stationnames={stations.station};
   plh0=cell2mat({stations.plh0}');
end

%xyz=plh2xyz([ plh0(:,1:2)*pi/180 plh0(:,3)]);
xyz=plh2xyz(plh0);
rdnap=etrs2rdnap(xyz);

%xywindow=[ 230 270 560 610]
xywindow=[ floor(min(rdnap(:,1))/5000)*5  ceil(max(rdnap(:,1))/5000)*5 ...
           floor(min(rdnap(:,2))/5000)*5  ceil(max(rdnap(:,2))/5000)*5 ] ...
           + [-10 5 -10 5];
xylegend=xywindow([ 1 3])+[ 5 5];

fprintf('\n      Lat [deg]  Lon [deg]  hght [m]\n')
for i=1:nstations
  fprintf('%s  %9.6f  %9.6f %8.3f\n',stationnames{i},[ plh0(i,1:2)*180/pi plh0(i,3) ])
end
fprintf('\n         X-RD [m]     Y-RD [m]   NAP [m]\n')
for i=1:nstations
  fprintf('%s  %12.3f %12.3f  %7.3f\n',stationnames{i},rdnap(i,:))
end

% Get screensize and position figure [ left bottom width height ]

scrsz = get(0,'ScreenSize');
%figpos=[ scrsz(1)+scrsz(3)*0.05 scrsz(2)+scrsz(4)*0.05 scrsz(3)*0.9 scrsz(4)*0.9]; 
figpos=[ scrsz(3)*0.05 scrsz(4)*0.05 scrsz(3)*0.05+scrsz(4)*0.94*8/10 scrsz(4)*0.94];

% Plot map

hf=figure('Name',[ project '_map_' item '_' id ],'NumberTitle','off','OuterPosition',figpos);
plot(rdnap(:,1)./1000,rdnap(:,2)./1000,'s');
hold on;
for i=1:nstations
  text(rdnap(i,1)/1000,rdnap(i,2)/1000+1.5,stationnames{i});
end


% Add items to the map

switch item
  case {'Steps','steps'}
    titlestr=[titlestr ' (Steps)'];
    for i=1:nstations
      if iscell(stations)
         station=stations{i};
         load([ station id '.mat']);
      else
         tseries=stations(i);
         station=tseries.station;
      end
      eventtype={ tseries.events.type }';
      eventdate=cell2mat({ tseries.events.year }');
      ll=0;
      for l=tseries.parindex(3)+1:tseries.parindex(4)
        ll=ll+1;
        tjump=tseries.tjump(ll);  
        idxevent=find(abs(eventdate-tjump) < 1/365);
        %if ~strcmpi(events{idxevent,2},'REF') || tjump < t1 || tjump > t2, continue; end
        if ~strcmpi(eventtype{idxevent},'REF'), continue; end
        step=[ tseries.xlat(l) tseries.xlon(l) tseries.xrad(l)  ];
        sigstep=sqrt([ tseries.qxlat(l,l) tseries.qxlon(l,l) tseries.qxrad(l,l)].*tseries.omt); 
        fprintf('%s  %8.3f  %4s  %9.2f %9.2f %9.2f  %6.2f %6.2f %6.2f\n',station,tjump,eventtype{idxevent},step*1000,sigstep*1000);
        pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[step(2) step(1)]*10000,'r','LineWidth',2);
        pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[0 step(3)]*5000,'k','LineWidth',2);
      end     
    end
    pltarrow(xylegend,xylegend+[.001 0]*10000,'r','LineWidth',2);
    pltarrow(xylegend,xylegend+[0 .002]*5000,'k','LineWidth',2);
    text(xylegend(1)+5,xylegend(2)+1,' 1 mm/y Horz','FontWeight','bold');
    text(xylegend(1),xylegend(2)+6,'2 mm/y Up','FontWeight','bold');
  case {'Velocity','velocity'}
    titlestr=[titlestr ' (Velocity)'];
    for i=1:nstations
      if iscell(stations)
         station=stations{i};
         load([ station id '.mat']);
      else
         tseries=stations(i);
         station=tseries.station;
      end
      vel=tseries.vel;
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[vel(2) vel(1)]*1000,'r','LineWidth',2);
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[0 vel(3)]*1000,'k','LineWidth',2);
    end
    pltarrow(xylegend,xylegend+[.005 0]*1000,'r','LineWidth',2);
    pltarrow(xylegend,xylegend+[0 .005]*1000,'k','LineWidth',2);
    text(xylegend(1)+5,xylegend(2)+1,' 5 mm/y Horz','FontWeight','bold');
    text(xylegend(1),xylegend(2)+6,'5 mm/y Up','FontWeight','bold');
  case {'Tempi','tempi'}
    titlestr=[titlestr ' (Temp.Influence)'];
    for i=1:nstations
      if iscell(stations)
         station=stations{i};
         load([ station id '.mat']);
      else
         tseries=stations(i);
         station=tseries.station;
      end
      ii=find(ismember(tseries.parlabel,'TempI'));
      if isempty(ii), continue; end
      tempi=[tseries.xlat(ii) tseries.xlon(ii) tseries.xrad(ii) ];
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[tempi(2) tempi(1)]*1667,'r','LineWidth',2);
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[0 tempi(3)]*1667,'k','LineWidth',2);
    end
    pltarrow(xylegend,xylegend+[.003 0]*1667,'r','LineWidth',2);
    pltarrow(xylegend,xylegend+[0 .003]*1667,'k','LineWidth',2);
    text(xylegend(1)+5,xylegend(2)+1,' 3 mm/10^0C Horz','FontWeight','bold');
    text(xylegend(1),xylegend(2)+6,'3 mm/10^0C Up','FontWeight','bold');
  case {'Annual','annual'}
    titlestr=[titlestr ' (Annual Term)'];
    for i=1:nstations
      if iscell(stations)
         station=stations{i};
         load([ station id '.mat']);
      else
         tseries=stations(i);
         station=tseries.station;
      end
      if tseries.parindex(2) == tseries.parindex(3), continue; end
      ii=tseries.parindex(2);
      phi=[9:3:360]'*pi/180;
      ellips=[ tseries.xlat(ii+1)*sin(phi)+tseries.xlat(ii+2)*cos(phi) ...
               tseries.xlon(ii+1)*sin(phi)+tseries.xlon(ii+2)*cos(phi) ...
               tseries.xrad(ii+1)*sin(phi)+tseries.xrad(ii+2)*cos(phi) ];
      amp3=sqrt(tseries.xrad(ii+1).^2+tseries.xrad(ii+2)^2);
      plot(rdnap(i,1)/1000+ellips(:,2)*1667,rdnap(i,2)/1000+ellips(:,1)*1667,'r','LineWidth',2);
      pltarrow(rdnap(i,1:2)/1000+ellips(end-1,[2 1])*1667,rdnap(i,1:2)/1000+ellips(end,[2 1])*1667,'r');
      plot(rdnap(i,1)/1000+[0;0],rdnap(i,2)/1000+[amp3 ;-amp3]*1667,'k','LineWidth',2);
      pltarrow(rdnap(i,1:2)/1000+[0 ellips(end-1,3)]*1667,rdnap(i,1:2)/1000+[0 ellips(end,3)]*1667,'k');
    end
    pltarrow(xylegend+[0 2],xylegend+[0 2]+[.003 0]*1667,'r','LineWidth',2);
    pltarrow(xylegend+[0 4],xylegend+[0 4]+[.003 0]*1667,'k','LineWidth',2);
    text(xylegend(1)+5,xylegend(2)+2,' 3 mm Horz','FontWeight','bold');
    text(xylegend(1)+5,xylegend(2)+4,' 3 mm Vert','FontWeight','bold');
  case {'AnnualP','annualp'}
    titlestr=[titlestr ' (Annual Phasor)'];
    for i=1:nstations
      if iscell(stations)
         station=stations{i};
         load([ station id '.mat']);
      else
         tseries=stations(i);
         station=tseries.station;
      end
      if tseries.parindex(2) == tseries.parindex(3), continue; end
      ii=tseries.parindex(2);
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[tseries.xlat(ii+2) tseries.xlat(ii+1)]*1000,'b','LineWidth',2);
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[tseries.xlon(ii+2) tseries.xlon(ii+1)]*1000,'r','LineWidth',2);
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[tseries.xrad(ii+2) tseries.xrad(ii+1)]*1000,'k','LineWidth',2);
    end
    pltarrow(xylegend,xylegend+[.005 0]*1000,'b','LineWidth',2);
    pltarrow(xylegend+[0 2],xylegend+[0 2]+[.005 0]*1000,'r','LineWidth',2);
    pltarrow(xylegend+[0 4],xylegend+[0 4]+[.005 0]*1000,'k','LineWidth',2);
    text(xylegend(1)+5,xylegend(2),' 5 mm Lan','FontWeight','bold');
    text(xylegend(1)+5,xylegend(2)+2,' 5 mm Lon','FontWeight','bold');
    text(xylegend(1)+5,xylegend(2)+4,' 5 mm Vert','FontWeight','bold');
  case {'AnnualA','annuala'}
    titlestr=[titlestr ' (Annual Amplitude)'];
    for i=1:nstations
      if iscell(stations)
         station=stations{i};
         load([ station id '.mat']);
      else
         tseries=stations(i);
         station=tseries.station;
      end
      if tseries.parindex(2) == tseries.parindex(3), continue; end
      amplitude=tseries.amplitude(:,1);
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[amplitude(2) amplitude(1)]*1000,'r','LineWidth',2);
      pltarrow(rdnap(i,1:2)/1000,rdnap(i,1:2)/1000+[0 amplitude(3)]*1000,'k','LineWidth',2);
    end
    pltarrow(xylegend,xylegend+[.005 0]*1000,'r','LineWidth',2);
    pltarrow(xylegend,xylegend+[0 .005]*1000,'k','LineWidth',2);
    text(xylegend(1)+5,xylegend(2)+1,' 5 mm Horz','FontWeight','bold');
    text(xylegend(1),xylegend(2)+6,'5 mm Up','FontWeight','bold');
  case {'Cov','cov'}
    titlestr=[titlestr ' (Emperical co-variance)'];
    for i=1:nstations
      if iscell(stations)
         station=stations{i};
         load([ station id '.mat']);
      else
         tseries=stations(i);
         station=tseries.station;
      end
      Q=cov(tseries.neuresidual);
      Q22=Q([2 1],[2 1]);
      sig3=sqrt(Q(3,3));
      pltellips(rdnap(i,1:2)/1000,Q22*1667.^2,'r','LineWidth',2);
      plot(rdnap(i,1)/1000+[0;0],rdnap(i,2)/1000+[sig3 ;-sig3]*1667,'k','LineWidth',2);
    end
    pltarrow(xylegend,xylegend+[.003 0]*1667,'r','LineWidth',2);
    pltarrow(xylegend,xylegend+[0 .003]*1667,'k','LineWidth',2);
    text(xylegend(1)+5,xylegend(2)+1,' 3 mm Horz','FontWeight','bold');
    text(xylegend(1),xylegend(2)+6,'3 mm Up','FontWeight','bold');
  case 'maponly'
  otherwise
    disp(['Invalid item ' item ]);
end
if ~isempty(id) && ~strcmp(item,'maponly')
  titlestr=[ titlestr(1:end-1) ', id=' id ')'];
end
title(titlestr,'Interpreter','none');
%set(gca,'YTickLabelRotation',90)
axis equal;
axis(xywindow);
xlabel('x RD [km]');
ylabel('Y RD [km]');
grid

if saveplots
  plotdir='plots';
  print(hf,'-dpng','-r600',fullfile(plotdir,[ 'map_' item id ' .png']))
end

end