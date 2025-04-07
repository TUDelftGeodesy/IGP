function [gpssecond,gpsweek]=date2gps(varargin)
%DATE2GPS    Convert Matlab date to gpssecond and gpsweek.
%   [GPSSECOND,GPSWEEK]=DATE2GPS(...) returns GPSSECOND with the GPS second 
%   in the GPS week GPSWEEK. Input is anything that is accepted by the 
%   matlab function datenum: date number, date string or date vector.
%   The returned GPSWEEK is always scalar, GPSSECOND is scalar or vector. 
% 
%   See also gps2date, date2mjd and datenum

%   (c) Hans van der Marel, Delft University of Technology, 2017.

%   Created:     7 August 2017 by Hans van der Marel
%   Modified:   

% Number of days since starting epoch of GPS weeks (Sunday 06-Jan-1980)

day0=datenum([1980 1 6 0 0 0]);
days = datenum(varargin{:}) - day0;
   
% Week number and second in week

gpsweek  = min(floor(days/7));
gpssecond = (days - gpsweek*7)*86400;

end