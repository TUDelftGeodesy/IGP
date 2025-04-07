function mdate=mjd2date(mjd)
%MJD2DATE    Convert Modified Julian Date to Matlab date number.
%   MDATE=MJD2DATE(MJD) returns the Matlab date number. Input is Modified 
%   Julian Date MJD (JD-2400000.5).
% 
%   See also date2mjd.

%   (c) Hans van der Marel, Delft University of Technology, 2016.

mdate=mjd + 678942;

end