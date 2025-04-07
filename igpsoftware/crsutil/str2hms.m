function hms=str2hms(str)
%STR2HMS Read hour, minute and second from string.
%   HMS=STR2HMS(STR) returns a three-element row vector HMS containing the 
%   decimal [ hour minute second]. Input is the time in string format 
%   'hh:mm:ss.sss'. The the hour can be outside the usual range of 0-23.
%
%   See also HMS2STR, STR2SEC, SEC2STR, HMS2STR and STR2HMS.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995-2013.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements
%                - Complete rewrite (note:cannot use datevec because
%                  this standard function is not very accurate for
%                  fractions of a second)

idx=findstr(str(1,:),':');
if length(idx) == 2
   hms=[str2num(str(:,1:idx(1)-1)) ...
     str2num(str(:,idx(1)+1:idx(2)-1))  ...
     str2num(str(:,idx(2)+1:end)) ];
elseif length(idx) == 1
   hms=[str2num(str(:,1:idx(1)-1)) ...
     str2num(str(:,idx(1)+1:end)) ...
     zeros(size(str,1),1)];
elseif length(idx) ==0 
   hms=sec2hms(str2num(str));
else
   hms=nan(size(str,3));
end
