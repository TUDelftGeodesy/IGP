function olc=olcencode(lat,lon,varargin)
%OLCENCODE   Encode Open Location Code (plus code).
%   OLC=OLCENCODE(LAT,LON) computes the ten most significant digits of
%   the Open Location Code OLC for points with latitude LAT and longitude
%   LON in degrees. LAT and LON can be scalar or arrays of equal length.
%   OLC is an character array with the 11 digit OLC codes, with a '+' as
%   separator between the 8th and 9th digit.
%
%   OLC=OLCENCODE(LAT,LON,LSD) adds LSD least significant digits to the
%   the ten most significant digits of the OLC code. LSD must be between
%   zero and five (default is 0).
%
%   The ten most significant digits represent an area of about 13.9 by 13.9
%   meter wide at the equator. Each LSD divides the area into a 4x5 grid 
%   and a single digit used to identify the grid square. A single grid 
%   refinement step reduces the area to approximately 3.5x2.8 meters in 
%   longitude and latitude direction; two to approximately 87x56 cm, three
%   to 22x11 cm, four to 5x2 cm and five to 14x4 mm. All at the equator.
%
%   If there are duplicate location codes, the function changes the '+'
%   separator to 'a', 'b', etc. A maximum of 26 duplicates can be 
%   handled.
%
%   (c) Hans van der Marel, Delft University of Technology, 2020

% Created:  26 August 2020 by Hans van der Marel
% Modified:

% Check the input arguments

if nargin < 2
    error('This function expects two inputs.')
end
if nargin ==3
    lsd=varargin{1};
    if lsd < 0 || lsd > 5 
        error('LSD must be between 0 and 5.')
    end
else 
    lsd=0;
end

% Symbol table

symb=['2' '3' '4' '5' '6' '7' '8' '9' 'C' 'F' 'G' 'H' 'J' 'M' 'P' 'Q' 'R' 'V' 'W' 'X' ]';

% Most significant 10 digits
%
%   The following provides an algorithm to encode the values from least 
%   significant digit to most significant digit:
% 
%   1.  Add 90 to the latitude and add 180 to the longitude, multiply both 
%       by 8000 and take the integer parts as latitude and longitude 
%       respectively
%   2.  Prefix the existing code with the symbol that has the integer part of 
%       longitude modulus 20
%   3.  Prefix the existing code with the symbol that has the integer part of 
%       latitude modulus 20
%   4.  Divide both longitude and latitude by 20
%   5.  Repeat from step 2 four more times.


plat=floor((lat(:)+90)*8000);
plon=floor(mod(lon(:)+180,360)*8000);

olc = [ symb(floor(mod(plat,20))+1) symb(floor(mod(plon,20))+1) ];
for k=1:4
   plat=plat/20;
   plon=plon/20;
   olc = [ symb(floor(mod(plat,20))+1) symb(floor(mod(plon,20))+1) olc ];
end

olc= [ olc(:,1:8) repmat('+',size(plat)) olc(:,9:10)];

if lsd > 0

% Least significant 5 digits
%
%   This differs from the above method in that each step produces a single 
%   character. This encodes latitude into base five and longitude into base four, 
%   and then multiplies the digits for each position together.
% 
%   The following provides an algorithm to encode the values from least 
%   significant digit to most significant digit:
% 
%   1. Add 90 to the latitude, multiply the fractional part by 2.5e7 and take 
%      the integer part as latitude
%   2. Add 180 to the longitude, multiply the fractional part by 8.192e6 and 
%      take the integer part as longitude
%   3. Multiply the integer part of longitude modulus 4 by the integer part 
%      of latitude modulus 5. NOTE: THE DESCRIPTION OF THIS STEP ON THE OFFICIAL
%      WEBSITE IS NOT CORRECT, IT SHOULD BE (CHECKED WITH THE JAVASCRIPT
%      CODE):
%      Multiply the integer part of latitude modulus 5 by a factor 4 and
%      add the integer part of longitude modulus 4.
%   4. Prefix the existing code with the symbol with the above value
%   5. Divide longitude by four and latitude by five
%   6. Repeat from step 2 four more times.

flat=floor(mod(lat(:)+90,1)*2.5e7);
flon=floor(mod(lon(:)+180,1)*8.192e6);

olclsd = symb( floor(mod(flon,4)) + floor(mod(flat,5))*4 + 1 );
for k=1:4
   flat=flat/5;
   flon=flon/4;
   olclsd = [ symb( floor(mod(flon,4)) + floor(mod(flat,5))*4 + 1 ) olclsd ];
end

olc = [ olc olclsd(:,1:lsd) ];

end

% Find duplicate olc and modify duplicates

for k=0:25
  [~, ind] = unique(olc, 'rows');
  duplicate_ind = setdiff(1:size(olc, 1), ind);
  if isempty(duplicate_ind), break; end
  % duplicate_value = olc(duplicate_ind, :);
  olc(duplicate_ind,9)= repmat(char(double('a')+k),size(duplicate_ind));
end

end
