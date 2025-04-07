function mdate=dyear2date(dyear)
%DYEAR2DATE    Convert decimal year to Matlab date number.
%   MDATE=DYEAR2DATE(DYEAR) returns the Matlab date number. Input is 
%   decimal year.
% 
%   See also date2year, date2mjd and mjd2date.

%   (c) Hans van der Marel, Delft University of Technology, 2016.


dmjd=51544.5+(dyear-2000)*365.25;
mdate=dmjd + 678942;

end