function str=hms2str(hms,fmt)
%HMS2STR Time to string conversion.
%   STR=HMS2STR(HMS) converts the row-vector HMS with decimal hour, minute
%   and second [ hour minute second] into the string format 'hh:mm:ss.sss' (
%   this is the default). 
%
%   STR=HMS2STR(HMS,FMT) replaces the default format by the optional format 
%   string FMT. The default format string is FMT='%02d:%02d:%06.3f'.  
%
%   The minute and second are converted to the usual range, the hour
%   can be outside the usual range of 0-23.
%
%   See also STR2HMS, HMS2SEC, SEC2HMS, SEC2STR and STR2SEC.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements
%                - Rewrite using num2str

if nargin==1, fmt='%02d:%02d:%06.3f';, end 

str=num2str(hms,fmt);