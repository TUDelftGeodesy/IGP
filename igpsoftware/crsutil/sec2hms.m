function hms=sec2hms(sec)
%SEC2HMS Time in seconds since 0hrs to hour, minutes and seconds.
%   HMS=SEC2HMS(SEC) returns a three-element row vector containing the 
%   decimal [ hour minute second]. Input is the time in seconds since midnight 
%   (0hrs). If the seconds are not within the usual range of one day the 
%   function returns a hour outside the range of 0-23.
%
%   See also HMS2SEC, HMS2STR, SEC2STR, STR2HMS and STR2SEC.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements
%               31 March 2015 by Hans van der Marel
%                - use fix instead of floor in computations

hour    = fix(sec/3600);
minute  = fix(sec/60);
minute  = minute - hour*60;
second  = sec - minute*60 - hour*3600;
hms = [hour minute second];

