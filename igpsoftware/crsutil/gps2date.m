function mdate=gps2date(gpssecond,gpsweek)
%GPS2DATE    Convert GPS second of week and GPS week to Matlab date number.
%   MDATE=GPS2DATE(GPSSECOND,GPSWEEK) returns the Matlab date number. Input 
%   is a scalar or vector GPSSECOND with the GPS second in the GPS week GPSWEEK.
% 
%   See also date2gps and mjd2date.

%   (c) Hans van der Marel, Delft University of Technology, 2017.

%   Created:     7 August 2017 by Hans van der Marel
%   Modified:   

day0=datenum([1980 1 6 0 0 0]);

mdate=day0+double(gpsweek)*7+gpssecond/86400;


end