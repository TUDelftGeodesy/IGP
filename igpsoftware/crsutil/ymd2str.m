function str=ymd2str(ymd)
%YMD2STR Date to string conversion.
%   STR=YMD2STR(YMD) converts the row-vector with decimal [year month day]
%   into the string format 'dd-mmm-yy'.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995-2013.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements
%                - Rewrite using datestr and datenum

str=datestr(datenum(ymd),'dd-mmm-yyyy');
