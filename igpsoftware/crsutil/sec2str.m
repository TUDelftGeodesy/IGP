function str=sec2str(sec,fmt)
%SEC2STR Time to string conversion.
%   STR=SEC2STR(HMS) converts the time in seconds since 0hrs into the 
%   string format 'hh:mm:ss.sss' (this is the default). The 
%   default format can be changed by including the optional format 
%   specifier SEC2STR(SEC,FMT). Default is FMT='%2d:%2d:%f6.3'.  
%   The resulting hour can be outside the usual range of 0-23.
%
%   See also STR2SEC, STR2HMS, HMS2STR, SEC2HMS and HMS2SEC.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements

hms=sec2hms(sec);
if nargin==1
  str=hms2str(hms);
else
  str=hms2str(hms,fmt);
end      
