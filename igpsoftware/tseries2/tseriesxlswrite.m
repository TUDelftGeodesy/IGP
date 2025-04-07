function tseriesxlswrite(stations,item,id,xlsfile)
%TSERIESXLSWRITE  Write timeseries to an excel file.
%   TSERIESXLSWRITE(TSERIES,ITEM,ID,XLSFILE) writes the timeseries components
%   in the cell array ITEM for the stations in the structure array TSERIES
%   to an excel file. ID is an optional identifier
%   and XLSFILE an optional excel filename.
%
%   TSERIESXLSWRITE(STATIONS,ITEM,ID,XLSFILE) writes the timeseries components
%   in the cell array ITEM for the stations in the cell array STATIONS
%   to an excel file. ID is the optional identifier 
%   of the fit (default '_fit'). The timeseries data must be stored in an 
%   mat-file STATION_ID.mat (_ is part of the ID). XLSFILE is an optional
%   output file name.
%
%   Default for ITEM is {'trend' 'residual'}, default for ID is '_fit' and
%   default for XLSFILE is "neu_<ITEMLIST>_<ID>.xlst".
%
%   Examples:
%      tseriesxlswrite(tseries);
%      tseriesxlswrite(stations);
%      tseriesxlswrite(stations,{'trend' 'residual'},'_fit');
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016.

% Check arguments

if nargin < 1, error('too few arguments (specify at least stationames');end
if ~iscell(stations) && ~isstruct(stations)
   error('first argument must be a cell array with station names or a structure array of timeseries.') 
end
if nargin < 2
   item={'trend' 'residual'};
end
item=cellstr(item);
if nargin < 3
   id='_fit';
end
if nargin < 4
   itemlist=item{1};
   for k=2:length(item)
     itemlist=[ itemlist '+' item{k} ];
   end
   xlsfile=[ 'neu_' itemlist id '.xlsx' ];
end

% Prepare

remjump=false;
if ismember(item,{'raw' 'fit'})
  remjump=true;
end

for k=1:length(item)
  if strcmp(item{k},'raw')
    item{k}='neu';
  else
    item{k}= [ 'neu' item{k} ];
  end
end

% Write all stations to excel file, one sheet per station

nstations=length(stations);    

header={'date' 'dyear' 'dN [mm]' 'dE [mm]' 'dU [mm]'};

for i=1:nstations

  if iscell(stations)
     station=stations{i};
     load([ station id '.mat']);
  elseif isstruct(stations)
    tseries=stations(i);
    station=tseries.station;
  end
  
  year=tseries.year;
  data=tseries.(item{1});
  if remjump, data=data-tseries.neujumps; end
  for k=2:length(item)
     data=data+tseries.(item{k});
  end          

  YMD=dyear2ymd(year);
  sdate=cellstr(datestr(YMD,'yyyy-mm-dd HH:MM'));
  
  result=[ year data.*1000];
  xlswrite(xlsfile,header,station,'A1')
  xlswrite(xlsfile,sdate,station,'A2')
  xlswrite(xlsfile,result,station,'B2')

end

end