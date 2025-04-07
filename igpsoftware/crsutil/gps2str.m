function str=gps2str(gpssecond,gpsweek);
%GPS2STR Convert gpsweek/seconds into a string
%   STR=GPS2STR(GPSSECOND,GPSWEEK) converts the vector GPSSECOND into a string 
%   of the form 'dd-mmm-yyyy HH:MM:SS.SSS'. GPSWEEK is the GPS weeknumber.
%
%   (c) Peter Joosten, Hans van der Marel, Delft University of Technology, 2000.

%   Created:    13 March 2000 by Peter Joosten and Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements
%                - Removed global variable GPSWEEK
%                - Complete rewrite using standard Matlab datestr function

day0=datenum([1980 1 6 0 0 0]);
day=floor(gpssecond/86400);
sec=gpssecond-day*86400;
day=day0+gpsweek*7+day;
imin=floor(sec/60);
sec=sec-imin*60;

str=[datestr(day+imin/1440,'dd-mmm-yyyy HH:MM') num2str(sec,':%06.3f') ];

