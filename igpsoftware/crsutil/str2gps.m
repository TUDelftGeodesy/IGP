function [gpssecond,gpsweek] = str2gps(str);
%STR2GPS Convert date/time string to GPS Week and Second past Sunday 0hrs. 
%   [GPSSECOND,GPSWEEK]=STR2GPS(STR) converts a time and date, given as string 
%   in the form of 'dd-mon-yyyy hh:mm:ss.sss' into the GPSWEEK and GPSSECOND 
%   into the GPSWEEK. GPSSECOND can be a vector, but GPSWEEK will always be 
%   scalar. 
%
%   (c) Hans van der Marel, Delft University of Technology, 2000.

%   Created:    13 March 2000 by Hans van der Marel and Peter Joosten
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements
%                - Removed global variable GPSWEEK
%                - Complete rewrite (note:cannot use datevec because
%                  this standard function is not very accurate for
%                  fractions of a second)

day0=datenum([1980 1 6 0 0 0]);

% split the string in a date and time part

idx=findstr(str(1,:),' ');
str1=str(:,1:idx(1)-1);
str2=str(:,idx(1)+1:end);

% convert both parts

gpsday=round(datenum(str1)-day0);
gpssecond=str2sec(str2);

gpsweek=min(floor(gpsday/7));
gpssecond = gpssecond + (gpsday - gpsweek*7) * 86400;

