function sec=str2sec(str)
%STR2SEC Convert time in string-format to seconds since 0hrs.
%   SEC=STR2SEC(STR) returns the time in seconds since midnight (0hrs). 
%   Input is the time in string format 'hh:mm:ss.sss'. The the hour
%   can be outside the usual range of 0-23.
%
%   See also STR2HMS, HMS2SEC, SEC2STR, SEC2HMS and HMS2STR.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements

hms=str2hms(str);
sec=hms2sec(hms);
