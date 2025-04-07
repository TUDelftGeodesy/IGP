function ymd=dyear2ymd(dyear)
%DYEAR2YMD Decimal Year to Gregorian Calendar.
%   YMD=DYEAR2YMD(DYEAR) returns the six-element row vector YMD with the 
%   [ year month day hour minute second] in the Gregorian Calendar. Input 
%   is the Decimal Year. 
%
%   See also YMD2DYEAR and MJD2YMD and SEC2HMS.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015.

%   Created:    31 March 2015 by Hans van der Marel
%   Modified:   

dmjd=51544.5+(dyear-2000)*365.25;
mjd=fix(dmjd+1e-7);
sec=(dmjd-mjd)*86400;

ymd= [ mjd2ymd(mjd) sec2hms(sec) ];


