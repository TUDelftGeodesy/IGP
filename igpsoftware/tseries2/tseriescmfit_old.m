function cm=tseriescmfit(stations,id)
%TSERIESCMFIT  Compute the common mode fit for GPS timeseries.
%   CM=TSERIESCMFIT(TSERIES,ID) computes the common mode CM for the timeseries 
%   in the structure array TSERIES. ID is an optional identifier. The
%   common mode is the weighted average (in least squares sense) of the 
%   temperature influence, atmospheric loading and harmonics of the 
%   indivudual series. The output CM is a structure with the estimated
%   common mode.
%
%   CM=TSERIESCMFIT(STATIONS,ID) computes the common mode CM of the timeseries 
%   for the stations in the cell array STATIONS.  ID is the optional 
%   identifier of the fit (default '_fit'). The timeseries data must be 
%   stored in an mat-file STATION_ID.mat (_ is part of the ID). 
%
%   Examples:
%      cm=tseriescmfit(tseries,'_fit');
%      cm=tseriescmfit(stations,'_fit');
%
%   See also TSERIESCMEVAL.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016.

% Check arguments

if nargin < 1, error('too few arguments (specify at least stationames');end
if ~iscell(stations) && ~isstruct(stations)
    error('first argument must be a cell array with station names or a structure array of timeseries.') 
end
nstations=length(stations);          
if nargin < 2
  id='_fit';
end

% Compute common mode for the temperature influence, atmospheric loading
% and harmonics

init=true;
env0=zeros(1,2);
for i=1:nstations
   if iscell(stations)
     station=stations{i};
     load([ station id '.mat']);
   elseif isstruct(stations)
     tseries=stations(i);
     station=tseries.station;
   end
   stationlist{i}=station;
   env0=env0+tseries.env0;
   % temperature influence, atmospheric loading and harmonics
   i1=tseries.parindex(1);
   i2=tseries.parindex(2);
   i3=tseries.parindex(3);
   if init
      harmonic=tseries.harmonic;
      parlabel=tseries.parlabel(i1+1:i3);
      parameters='';
      for j=i1+1:i3
         parameters=[parameters sprintf(' %7s',tseries.parlabel{j}) ];
      end
      xlat=zeros(i3-i1,1);
      xlon=zeros(i3-i1,1);
      xrad=zeros(i3-i1,1);
      xlatstations=zeros(i3-i1,nstations);
      xlonstations=zeros(i3-i1,nstations);
      xradstations=zeros(i3-i1,nstations);     
      qxlat=zeros(i3-i1,i3-i1);
      qxlon=zeros(i3-i1,i3-i1);
      qxrad=zeros(i3-i1,i3-i1);
      init=false;
   else
      if i3-i1 ~= size(xlat,1), error('all series must have similar parameters'); end
   end
   xlati=tseries.xlat(i1+1:i3);
   xloni=tseries.xlon(i1+1:i3);
   xradi=tseries.xrad(i1+1:i3);
   qxlati=tseries.qxlat(i1+1:i3,i1+1:i3); 
   qxloni=tseries.qxlon(i1+1:i3,i1+1:i3);
   qxradi=tseries.qxrad(i1+1:i3,i1+1:i3);
   wlati=inv(qxlati);
   xlat=xlat+wlati*xlati;
   qxlat=qxlat+wlati;
   wloni=inv(qxloni);
   xlon=xlon+wloni*xloni;
   qxlon=qxlon+wloni;
   wradi=inv(qxradi);
   xrad=xrad+wradi*xradi;
   qxrad=qxrad+wradi;
   xlatstations(:,i)=xlati;
   xlonstations(:,i)=xloni;
   xradstations(:,i)=xradi;
end
qxlat=inv(qxlat);
xlat=qxlat*xlat;
qxlon=inv(qxlon);
xlon=qxlon*xlon;
qxrad=inv(qxrad);
xrad=qxrad*xrad;

env0=env0./nstations;

% Print common mode parameters

fmt2='%s  '; 
for j=1:size(xlatstations,1)
  fmt2=[fmt2 ' %7.2f'];
end
fmt2=[fmt2 '\n'];

fprintf('\nLat\n')
fprintf('       %s\n',parameters)
for i=1:nstations
   fprintf(fmt2,stationlist{i},xlatstations(:,i)*1000);
end
fprintf('\n')
fprintf(fmt2,'    ',xlat*1000);
%fprintf(fmt2,'    ',sqrt(diag(qxlat))*1000);

fprintf('\nLon\n')
fprintf('       %s\n',parameters)
for i=1:nstations
   fprintf(fmt2,stationlist{i},xlonstations(:,i)*1000);
end
fprintf('\n')
fprintf(fmt2,'    ',xlon*1000);
%fprintf(fmt2,'    ',sqrt(diag(qxlon))*1000);

fprintf('\nRad\n')
fprintf('       %s\n',parameters)
for i=1:nstations
   fprintf(fmt2,stationlist{i},xradstations(:,i)*1000);
end
fprintf('\n')
fprintf(fmt2,'    ',xrad*1000);
%fprintf(fmt2,'    ',sqrt(diag(qxrad))*1000);

cm.stations=stationlist;
cm.parlabel=parlabel;
cm.harmonic=harmonic;

cm.env0=env0;

cm.xlat=xlat;
cm.xlon=xlon;
cm.xrad=xrad;
cm.qxlat=qxlat;
cm.qxlon=qxlon;
cm.qxrad=qxrad;
cm.xlatstations=xlatstations;
cm.xlonstations=xlonstations;
cm.xradstations=xradstations;

end