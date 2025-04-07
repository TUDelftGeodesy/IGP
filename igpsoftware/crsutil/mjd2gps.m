function [gpssecond,gpsweek]=mjd2gps(mjd)
%MJD2GPS Modified Julian Date (MJD) to GPS Week and Second. 
%   [GPSSECOND,GPSWEEK]=MJD2YMD(MJD) returns the GPSSECOND and GPSWEEK for
%   Modified Julian Date MJD (JD-2400000.5).GPSWEEK is the weeknumber since
%   Sunday, 6th January 1980, GPSSECOND is the number of seconds past midnight 
%   of last saterday/sunday referenced by GPSWEEK. GPSSECOND can be a vector, 
%   but GPSWEEK will always be scalar.
%
%   See also GPS2MJD, MJD2YMD and YMD2MJD.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Changed function output from 1 to 2 arguments
%                - Updated description and copyright statements

% Number of days since starting epoch of GPS weeks (Sunday 06-Jan-1980)
mjds = mjd - 44244;
      
% Week number and second in week

gpsweek  = min(floor(mjds/7));
gpssecond = (mjds - gpsweek*7)*86400;
