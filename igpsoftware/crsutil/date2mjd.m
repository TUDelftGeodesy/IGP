function mjd=date2mjd(varargin)
%DATE2MJD    Convert Matlab date to Modified Julian Date
%   MJD=DATE2MJD(...) returns the Modified Julian Date MJD (JD-2400000.5).
%   The input can be anything accepted by the matlab function datenum.
% 
%   See also datenum, ymd2mjd and mjd2date.

%   (c) Hans van der Marel, Delft University of Technology, 2016.

% mdate=datenum(varargin{:});
% v=datevec(mdate);
% mjd = ymd2mjd(v(:,1:3)) + ( v(:,4) + v(:,5)/60 + v(:,6)/3600 ) / 24;

mjd=datenum(varargin{:}) - 678942;

end