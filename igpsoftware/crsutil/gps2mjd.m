function mjd=gps2mjd(gpssecond,gpsweek)
%GPS2MJD  Convert GPS Week and Second to Modified Julian Date (MJD). 
%   MDJ=GPS2MJD(GPSSECOND,GPSWEEK) returns the Modified Julian Date MJD 
%   (JD-2400000.5) for GPSWEEK and GPSSECOND. GPSWEEK is the weeknumber since
%   Sunday, 6th January 1980, GPSSECOND is the number of seconds past midnight 
%   of last saterday/sunday.
%
%   See also MJD2GPS, MJD2YMD and YMD2MJD.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    39 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Changed function call from 1 to 2 arguments
%                - Updated description and copyright statements

mjd = gpsweek*7 + gpssecond/86400 + 44244;
