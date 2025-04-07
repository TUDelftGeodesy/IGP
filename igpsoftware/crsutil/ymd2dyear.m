function dyear=ymd2dyear(ymd)
%YMD2DYEAR Gregorian Calendar to Decimal Year
%   DYEAR=YMD2DYEAR(YMD) returns the Decimal Year at 0h. Input is a 
%   three-element row vector YMD containing the year, month and day 
%   [ year month day ] in the Gregorian Calendar, or a string in the 
%   format 'yy-mm-dd' or 'dd-mmm-yy'. Acceptable years are  
%     00-49, interpreted as 2000-2049,
%     50-99,     "       "  1950-1999,
%       100 upwards, interpreted literally.
%   If the day is not OK the function returns the MJD into the next
%   (previous) month or year. Fraction of days are carried into the DYEAR.
%   If the input year or month are not OK the function returns NaN. 
%
%   DYEAR=YMD2DYEAR(DATEVEC) returns the Decimal Year with as input a
%   6 element date vector with year, month, day, hour, minute and seconds.
%
%   See also DYEAR2YMD and YMD2MJD.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015.

%   Created:    31 March 2015 by Hans van der Marel
%   Modified:  

if ~isstr(ymd) && size(ymd,2) == 6
  fracday=(ymd(:,4)+ymd(:,5)./60+ymd(:,6)./3600)./24;
  ymd=ymd(:,1:3);
else
  fracday=0;
end

dyear=(ymd2mjd(ymd)-51544.5+fracday)./365.25+2000;

end