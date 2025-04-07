function [temperature,temperature_week] = get_meteo(dates,station,start,stop)

dates = datenum(dates,'yyyymmddHH');
Ndates = size(dates,1);

%download_file = ['knmi_download_stn' station '_' start '_' stop '.txt'];
download_file = 'knmi_download_stn280_2008010100_2019030321.txt';
call = ['wget -O ' download_file ' --post-data="stns=' station '&start=' start '&end=' stop '" http://projects.knmi.nl/klimatologie/uurgegevens/getdata_uur.cgi'];
%system(call);

meteo = textread(download_file,'','headerlines',34,'delimiter',',');

meteo_date = num2str(meteo(:,2));
meteo_hour = num2str(meteo(:,3),'%02d');

temperature = NaN(Ndates,1);
temperature_week = NaN(Ndates,1);
for v = 1:Ndates
  date = datestr(dates(v),'yyyymmdd');
  hour = datestr(dates(v),'HH');
    
  idx1 = strmatch(date,meteo_date);
  idx2 = strmatch(hour,meteo_hour(idx1,:));
  temperature(v) = 0.1*meteo(idx1(idx2),8);
  temperature_week(v) = 0.1*mean(meteo(idx1(idx2)-7*24:idx1(idx2),8));
end

