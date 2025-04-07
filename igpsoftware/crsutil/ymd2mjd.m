function mjd=ymd2mjd(ymd)
%YMD2MJD Gregorian Calendar to Modified Julian Date
%   MJD=YMD2MJD(YMD) returns the Modified Julian Date MJD (JD-2400000.5) 
%   for 0 hrs. Input is a three-element row vector YMD containing the 
%   decimal year, month and day [ year month day ] in the Gregorian Calendar,
%   or a string in the format 'yy-mm-dd' or 'dd-mmm-yy'.         
%   Acceptable years are  
%     00-49, interpreted as 2000-2049,
%     50-99,     "       "  1950-1999,
%       100 upwards, interpreted literally.
%   If the day is not OK the function returns the MJD into the next
%   (previous) month or year. Fraction of days are carried into the mjd.
%   If the input year or month are not OK the function returns NaN. 
%   Fractions of year and month are meaningless.
%
%   See also MJD2YMD, YMD2STR and STR2YMD.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Vectorized the default century part
%                - Updated description and copyright statements

if isstr(ymd)
  ymd2=str2ymd(ymd);
  year=round(ymd2(:,1));
  month=round(ymd2(:,2));
  day=ymd2(:,3);
else
  year=round(ymd(:,1));
  month=round(ymd(:,2));
  day=ymd(:,3);
end

% Add default century 
idx=(year>=0 & year<=49);
year(idx)=year(idx)+2000;
idx=(year>=50 & year<=99);
year(idx)=year(idx)+1900;

% Validate year and month
idx=find(year<-4699 | month<1 | month>12);
year(idx)=NaN;
month(idx)=NaN;

% Modified Julian Date
mjd=   fix((1461*(year-fix((12-month)/10)+4712))/4) ... 
   +   fix((306*rem(month+9,12)+5)/10) ...
   -   fix((3*fix((year-fix((12-month)/10)+4900)/100))/4) ... 
   +   day - 2399904;


