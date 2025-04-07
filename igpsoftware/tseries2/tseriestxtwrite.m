function tseriestxtwrite(stations,item,id,basename)
%TSERIESTXTWRITE  Write timeseries to an ascii text file.
%   TSERIESTXTWRITE(TSERIES,ITEM,ID) writes the timeseries components
%   in the cell array ITEM for the stations in the structure array TSERIES
%   to an ascii text file for each station. ID is an optional identifier.
%
%   TSERIESTXTWRITE(STATIONS,ITEM,ID) writes the timeseries components
%   in the cell array ITEM for the stations in the cell array STATIONS
%   to an ascii text file for each station. ID is the optional identifier 
%   of the fit (default '_fit'). The timeseries data must be stored in an 
%   mat-file STATION_ID.mat (_ is part of the ID).
%
%   The output files are called neu_<STATION>_<ITEMLIST>_<ID>.txt
%
%   Examples:
%      tseriestxtwrite(tseries,'raw');
%      tseriestxtwrite(stations,{'harmonic' 'tempi'},'_fit');
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
   basename= 'neu';
end

% Prepare

remjump=false;
if ismember(item,{'nostep' 'fit'})
  remjump=true;
end

itemlist=item{1};
for k=2:length(item)
  itemlist=[ itemlist '+' item{k} ];
end

for k=1:length(item)
  %if strcmp(item{k},{'raw' 'nostep'})
  if ismember(item{k},{'raw' 'nostep'})
    item{k}='neu';
  else
    item{k}= [ 'neu' item{k} ];
  end
end


% Write all stations to excel file, one sheet per station

nstations=length(stations);    

plh0=nan(nstations,3);
scor=zeros(nstations,6);

for i=1:nstations

  % Get the data
    
  if iscell(stations)
    station=stations{i};
    load([ station id '.mat']);
  elseif isstruct(stations)
    tseries=stations(i);
    station=tseries.station;
  end
  stationlist{i}=station;
  
  year=tseries.year;
  data=tseries.(item{1});
  if remjump, data=data-tseries.neujumps; end
  for k=2:length(item)
     data=data+tseries.(item{k});
  end         
  
  % Save to file
  
  YMD=dyear2ymd(year);
  sdate=cellstr(datestr(YMD,'yyyy mm dd HH MM'));
  
  fid=fopen([ basename '_' station '_' itemlist id '.txt'],'w');
  fprintf(fid,'year mm dd HH MM   dN [mm]   dE [mm]   dU [mm]\n');
  for k=1:numel(year)
    fprintf(fid,'%16s %9.2f %9.2f %9.2f\n',sdate{k},data(k,:)*1000);
  end
  fclose(fid);

  % Collect station coordinates and compute station co-variance

  if isfield(tseries,'plh0')
     plh0(i,:)=tseries.plh0;
  end
  Q=cov(tseries.neuresidual);
  scor(i,:)=covreformat(Q,'qmat','scor'); %  sx, sy, sz, rxy, rxz, ryz

end

% In previous versions plh0 was not saved, so may have to get is from
% somewhere else

if all(isnan(plh0(:)))
   for i=1:nstations
      station=stations{i};
      load([ station '.mat']);
      plh0(i,:)=tseries.plh0;
   end
end

% create file with coordinates and covariance

%xyz=plh2xyz([ plh0(:,1:2)*pi/180 plh0(:,3)]);
xyz=plh2xyz(plh0);
rdnap=etrs2rdnap(xyz,'zero');

fid=fopen(['station_crd' id '.txt'],'w');
fprintf(fid,'      Lat [deg]  Lon [deg]  hght [m]      X-RD [m]     Y-RD [m]   NAP [m]\n');
for i=1:nstations
  station=stationlist{i};
  %fprintf(fid,'%s  %9.6f  %9.6f %9.3f  %12.3f %12.3f  %8.3f\n',station,plh0(i,:),rdnap(i,:));
  fprintf(fid,'%s  %9.6f  %9.6f %9.3f  %12.3f %12.3f  %8.3f\n',station,plh0(i,1:2)*180/pi,plh0(i,3),rdnap(i,:));
end
fclose(fid);

scor(:,1:3)=scor(:,1:3)*1000;

fid=fopen(['station_cov' id '.txt'],'w');
fprintf(fid,'       sN [mm]  sE [mm]  sU [mm]  rNE [-]  rNU [-]  rEU [-]\n');
for i=1:nstations
  station=stationlist{i};
  fprintf(fid,'%s  %8.2f %8.2f %8.2f %8.2f %8.2f %8.2f\n',station,scor(i,:));
end
fclose(fid);

end

