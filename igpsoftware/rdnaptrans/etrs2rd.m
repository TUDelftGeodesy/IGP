function rd=etrs2rd(latlon,hmsl)
%ETRS2RD  Convert ETRS89 latitude and longitude into RD coordinates.
%  RD=ETRS2RD(LATLON) convert ETRS89 latitude and longitude in LATLON into
%  output matrix RD with Dutch RD coordinates. LATLON is a matrix with in
%  the first column the ETRS89 latitude in degrees, and in its second
%  column the longitude in degrees. RD has the same size as LATLON and 
%  contains the grid coordinates in meters. The height is assumed to be
%  zero in the Dutch NAP system.
%
%  RD=ETRS2RD(LATLON,HMSL) does the same, but uses HMSL for the height
%  in the Dutch NAP system. This height does not have to be very accurate,
%  and has only an effect at the mm level. HMSL can be scalar or an array.
%
%  Examples:
%  
%    rd=etrs2rd([51.8237   5.7931]);
%    rd=etrs2rd([50.763104 5.931024],135);
% 
%  See also rd2etrs, etrs2rdnap, rdnap2etrs, etrs2nap and nap2etrs.
% 
%  (c) Hans van der Marel, Delft University of Technology, 2004-2019

% Created:  30 Apr 2019 by Hans van der Marel, TUD
% Modified:  1 Oct 2023 by Hans van der Marel, TUD
%            - Set RD correction grid extension to ZERO to avoid NaN's

% Check input arguments

if nargin < 1 || nargin > 2, error('Function needs one or two input arguments.'); end
if nargin < 2 
    hmsl=0;
end

% Compute rd coordinates using approximate geoid height of 44 m

hell=hmsl+repmat(44,[size(latlon,1),1]);
rdnap1=etrs2rdnap([ latlon*pi/180 hell],'plh','ZERO');

% Compute rd coordinates again, but corrected now using actual geoid height

hell=hell-rdnap1(:,3)+hmsl;
rdnap2=etrs2rdnap([ latlon*pi/180 hell],'plh','ZERO');

rd=rdnap2(:,1:2);

end