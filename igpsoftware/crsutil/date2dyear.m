function dyear=date2dyear(varargin)
%DATE2DYEAR  Convert Matlab date to Decimal year
%   DYEAR=DATE2DYEAR(...) returns the decimal year. The input can be 
%   anything accepted by the matlab function datenum.
% 
%   See also datenum, dyear2date, date2mjd and mjd2date.

%   (c) Hans van der Marel, Delft University of Technology, 2016.

mjd=datenum(varargin{:}) - 678942;
dyear=(mjd-51544.5)./365.25+2000;

end