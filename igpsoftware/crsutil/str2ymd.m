function ymd=str2ymd(str)
%STR2YMD Convert date in string format into decimal form.
%   YMD=STR2YMD(STR) returns the three-element row vector with the decimal
%   [year month day ]. Input is a string STR with the date given as 'dd-mm-yy' 
%    or 'dd-mmm-yy', with the month numeric or as Jan, Feb, etc.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995-2013.

%   Created:    29 April 1995 by Hans van der Marel
%   Modified:   13 June 2013 by Hans van der Marel
%                - Updated description and copyright statements
%                - Rewrite using datevec

% Determine date format

idx=strfind(str(1,:),'-');
if length(idx) ~= 2
  error('input is not a valid date format')
end
if isletter(str(1,idx(1)+1:idx(2)-1))
  fmt='dd-mmm-yyyy';
elseif length(str(1,1:idx(1)-1)) == 4
  fmt='yyyy-mm-dd';
elseif length(str(1,idx(2)+1:end)) == 4
  fmt='dd-mm-yyyy';
else
  fmt='dd-mm-yy';
end

% Convert

[y m d]=datevec(str,fmt);
ymd=[ y m d];