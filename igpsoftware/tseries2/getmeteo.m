function [xday,Tday,Pday,xhour,Thour,Phour]=getmeteo(stationname,files,varargin)
%GETMETEO   Read temperature and pressure data from KNMI hourly meteo files.
%   GETMETEO(STATIONNAME,FILES) reads the temperature and pressure data
%   for station STATIONNAME from the KNMI "uurgegevens" files with file
%   names in the cell array FILES. The data is saved to a mat file with
%   the name meteo_<STATIONNAME>.mat. Several plots of the meteo data
%   are produced.
%
%   [XDAY,TDAY,PDAY,XHOUR,THOUR,PHOUR]=GETMETEO(STATIONNAME,FILES) outputs
%   the data as variables, no data is saved unless the dosave is set.
%
%   [...]=GETMETEO(STATIONNAME,FILES,OPTION,VALUE,...) accept the following
%   option value pairs
%
%     doplot   true|false   Plot the meteo data (default true)
%     dosave   true|false   Save the data to disk *)
%     saveplot true|false   Save plots to disk (default true)
%
%   *) default is true if no output is saved to variables, false otherwise
%
%   Example:
% 
%   unzip('http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/uurgeg_280_2001-2010.zip','meteo');
%   unzip('http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/uurgeg_280_2011-2020.zip','meteo');
%   getmeteo('Eelde',{'meteo/uurgeg_280_2001-2010.txt';'meteo/uurgeg_280_2011-2020.txt'}); 
%
%   (c) H. van der Marel, Delft University of Technology, 2015-2016, 2021.

% Check input arguments

if nargin < 2
   error('This function needs a station and file names as input')
end

% Default options

opt.doplot=true;
opt.saveplot=true;
opt.dosave=true;
if nargout > 1
  opt.dosave=false;    
end

% Check option/value pairs

for k=1:2:numel(varargin)-1
   if isfield(opt,varargin{k})
      opt.(varargin{k})=varargin{k+1};
   else
      error(['Illegal option ' varargin{k} ])
   end
end

% Read the data from file

%  http://www.knmi.nl/nederland-nu/klimatologie/uurgegevens

url='http://www.knmi.nl/nederland-nu/klimatologie/uurgegevens';

% files = { 'meteo/uurgeg_380_1991-2000.txt'  ; ... 
%           'meteo/uurgeg_380_2001-2010.txt'  ; ...
%           'meteo/uurgeg_380_2011-2020.txt'  };
% files = { 'meteo/uurgeg_280_2011-2020.txt' }; 

YYYYMMDD=[];
        
xday=[];
Tday=[];
Pday=[];
xhour=[];
Thour=[];
Phour=[];

for k=1:size(files,1)
   
   meteo=readmeteo(files{k});

   year=floor(meteo.YYYYMMDD ./ 10000);
   month=floor(meteo.YYYYMMDD ./ 100) -year*100;
   day=meteo.YYYYMMDD - year*10000 - month*100;

   xk=datenum(year,month,day,meteo.HH,0,0);
     
   xhour=[ xhour ; xk ];              % Matlab date number
   Thour=[ Thour ; 0.1*meteo.T ];     % Temperature [C]
   Phour=[ Phour ; 0.1*meteo.P ];     % Pressure [hPa]

   [Tdayk,xdaystr]=grpstats(0.1*meteo.T,meteo.YYYYMMDD,{'mean' 'gname'});
   [Pdayk]=grpstats(0.1*meteo.P,meteo.YYYYMMDD,{'mean'});
   xdayk=datenum(xdaystr,'yyyymmdd');

   YYYYMMDD=[ YYYYMMDD ; xdaystr ];
   
   xday=[xday; xdayk];
   Tday=[Tday; Tdayk];
   Pday=[Pday; Pdayk];

end

% Optional plots

if opt.doplot

  tmp=datevec(max(xday));lastyear=tmp(1);
  tmp=datevec(min(xday));firstyear=tmp(1);
    
  figure
  plot(xhour,Thour,'c')
  hold on
  plot(xday,Tday,'b')
  datetick('x')
  legend('Hourly Temp','Daily Average')
  ylabel('Temp [^oC]')
  title([ stationname ' KNMI'])
  a=axis();
  axis([a(1) datenum(lastyear+1,0,0) a(3:4) ])
  if opt.saveplot
    print('-dpng','-r 300',[ stationname '_Temp.png'] )
  end
  
  figure
  [year,~,~,]=datevec(xhour);
  plot(xhour-datenum(year,0,0),Thour,'c.')
  hold on
  [year,~,~,]=datevec(xday);
  iday=round(xday-datenum(year,0,0));
  plot(iday,Tday,'b.')
  Tmean=grpstats(Tday,iday);
  plot(Tmean,'g-','LineWidth',2)
  legend('Hourly Temp','Daily Average','Mean')
  datetick('x','ddmmm')
  ylabel('Temp [^oC]')
  title([ stationname ' KNMI (' num2str(firstyear) '-' num2str(lastyear) ')'])
  if opt.saveplot
    print('-dpng','-r 300',[ stationname '_Yearly_Temp.png'] )
  end

  figure
  plot(xday,Pday)
  datetick('x')
  ylabel('Pressure [hPa]')
  title([ stationname ' KNMI'])
  a=axis();
  axis([a(1) datenum(lastyear+1,0,0) a(3:4) ])
  if opt.saveplot
    print('-dpng','-r 300',[stationname '_Pressure.png'] )
  end

  figure
  [year,~,~,]=datevec(xhour);
  plot(xhour-datenum(year,0,0),Phour,'c.')
  hold on
  [year,~,~,]=datevec(xday);
  iday=round(xday-datenum(year,0,0));
  plot(iday,Pday,'b.')
  Tmean=grpstats(Pday,iday);
  plot(Tmean,'g-','LineWidth',2)
  legend('Hourly Pressure','Daily Average','Mean')
  datetick('x','ddmmm')
  ylabel('Pressure [hPa]')
  title([ stationname ' KNMI (' num2str(firstyear) '-' num2str(lastyear) ')'])
  if opt.saveplot
    print('-dpng','-r 300',[ stationname '_Yearly_Pressure.png'] )
  end

end

% Save optionally to mat file

if opt.dosave
  YYYYMMDD=cell2mat(YYYYMMDD);
  save([ 'meteo_' stationname '.mat'],'xday','Tday','Pday','YYYYMMDD');  
end

end


function meteo=readmeteo(knmifile)

fid=fopen(knmifile);
while ~feof(fid)
     line=fgetl(fid);
     if strncmp(line,'# STN,YYYY',10), break, end
end
fields=textscan(line(2:end),'%s','delimiter',',');
fields=fields{1};
nfields=size(fields,1);

format='';
for k=1:nfields
  format=[format '%f '];
end

line=fgetl(fid);
C = textscan(fid, format, 'delimiter', ',','EmptyValue', -Inf);
fclose(fid);

for k=1:nfields
   meteo.(fields{k})=C{k};
end

end