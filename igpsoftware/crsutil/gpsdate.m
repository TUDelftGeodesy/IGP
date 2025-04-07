function gpsdate(str)
%GPSDATE  Display current date, GPSWeek, Day of Week and Day of Year
%   GPSDATE(DATESTR) prints the GPS week number, day of week and
%   day of year for the given date in DATESTR.
%
%   GPSDATE does the same for the current date.
%
%   (c) Hans van der Marel, Delft University of Technology, 2000.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements

if nargin<1
  str=date;
end

for l=1:size(str,1)

  ymd=str2ymd(str(l,:));
  mjd=ymd2mjd(ymd);
  [gpssecond,gpsweek]=mjd2gps(mjd);
  doy=mjd-ymd2mjd([ ymd(1) 1 1])+1;

  disp([ datestr([ ymd [0 0 0]],1) ' -> week/dow ' num2str(gpsweek) '/' num2str(gpssecond/86400) ...
           ', doy ' sprintf('%3.0f',doy) ...
           ', gpsseconds ' sprintf('%8.1f',gpssecond) '-' sprintf('%8.1f',gpssecond+86400)])

end
