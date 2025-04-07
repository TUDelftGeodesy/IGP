function [s,f]=tseriesperiodogram(stations,item,id,stack)
%TSERIESPERIODOGRAM  Plot timeseries periodogram for one or more stations.
%   TSERIESPERIODOGRAM(TSERIES,ITEM,ID) plots the periodogram for the timeseries 
%   components in the cell array ITEM for the stations in the structure 
%   array TSERIES. ID is an optional identifier. The periodogram is
%   plotted.
%
%   TSERIESPERIODOGRAM(STATIONS,ITEM,ID) plots the periodogram for the timeseries 
%   components in the cell array ITEM for the stations in the cell array 
%   STATIONS.  ID is the optional identifier of the fit (default '_fit'). 
%   The timeseries data must be stored in an mat-file STATION_ID.mat 
%   (_ is part of the ID). The periodogram is plotted.
%
%   [S,F]=TSERIESPERIODOGRAM(...,STACK) is an exended calling sequence with
%   STACK an optional logical for stacking (Default true) and output of the 
%   periodogram [S,F]. If STACK is false multiple periodograms are plotted.
%
%   Examples:
%      tseriesperiodogram(tseries,'raw');
%      tseriesperiodogram(stations,{'harmonic' 'tempi'},'_fit');
%      tseriesplotcomponent(stations,{'harmonic' 'tempi'},'_fit','Ameland');
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016.

% Check arguments

if nargin < 2, error('too few arguments (specify at least stationames and item');end
if ischar(stations)
    stations=cellstr(stations);
end
if ~iscell(stations) && ~isstruct(stations)
    error('first argument must be a cell array with station names or a structure array of timeseries.') 
end
if nargin < 2
    item='raw';
end
item=cellstr(item);
if nargin < 3
    id='_fit';
end
if nargin < 4
   stack=true;
end

saveplots=false;

% Prepare

remjump=false;
if ismember(item,{'raw' 'fit'})
  remjump=true;
end

itemlist=[ item{1}];
for k=2:length(item)
  itemlist=[ itemlist ' + ' item{k} ];
end
if ~isempty(id)
    titlestr=[ '(' itemlist ', id=' id ')'];
else
    titlestr=[ '(' itemlist ')'];
end

for k=1:length(item)
  if strcmp(item{k},'raw')
    item{k}='neu';
  else
    item{k}= [ 'neu' item{k} ];
  end
end

% Compute periodogram

nstations=length(stations);          
for i=1:nstations
  if iscell(stations)
     station=stations{i};
     load([ station id '.mat']);
  elseif isstruct(stations)
     tseries=stations(i);
     station=tseries.station;
  end
  data=tseries.(item{1});
  if remjump, data=data-tseries.neujumps; end
  for k=2:length(item)
     data=data+tseries.(item{k});
  end
  if stack
    [s{i},f{i}]=plomb([data(:,1)*10 data(:,2)*100 data(:,3)*1000],tseries.year,[0.1:0.02:120]);
  else
    pltspectrum(tseries.year,data,[ 'Periodogram ' station ' ' titlestr ])
    if saveplots
       plotdir='plots';
       filename=['periodogram_' station '_' strrep(itemlist,' ','') id ]; 
       print(fullfile(plotdir,[ filename '.png']),'-dpng','-r600')
    end
  end
end

if stack
  stot=zeros(size(s{1}));
  for k=1:nstations
    stot=stot+s{k}; 
  end
  stot=stot./nstations;
  if nstations > 1
     pltspectrum(f{1},stot,[ 'Stacked periodogram ' titlestr],false)
     filename=['periodogram_stacked_' strrep(itemlist,' ','') id ]; 
  else
     pltspectrum(f{1},stot,[ 'Periodogram ' station ' ' titlestr],false)
     filename=['periodogram_' station '_' strrep(itemlist,' ','') id ]; 
  end
  if saveplots
     plotdir='plots';
     print(fullfile(plotdir,[ filename '.png']),'-dpng','-r600')
  end
end

end