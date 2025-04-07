function sec=hms2sec(hms)
%HMS2SEC Hour, minutes and seconds to time in seconds since 0hrs.
%   SEC=HMS2SEC(HMS) returns the time in seconds since midnight (0hrs). 
%   Input is a three-element row vector containing the decimal 
%   [ hour minute second] or a string in the format 'hh:mm:ss.sssss'. 
%   If the hour, minute or second are not within the usual range the function 
%   returns the corresponding time in seconds since 0hrs anyhow.
%
%   See also SEC2HMS, STR2HMS, HMS2STR.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements

if isstr(hms)
  hms=str2hms(hms);
end

sec=(hms(:,1)*60+hms(:,2))*60+hms(:,3);
