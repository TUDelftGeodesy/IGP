function [gmst,omegae]=ut2gmst(ut1,model)
%UT2GMST  Compute Greenwich Mean Siderial Time from UT1
%    GMST=UT2GMST(UT1) returns the Greenwich Mean Siderial Time GMST 
%    [0-2pi rad] for UT1. UT1 is a matlab date number or matlab date string
%    representing UT1.
%
%    [GMST,OMEGAE]=UT2GMST(UT1) returns also the rotation rate of the 
%    Earth [rev/day].
%
%    [...]=UT2GMST(...,MODEL) selects the computation method, possible choices
%    for MODEL are APPROXIMATE and IAU-82. The default is the IAU-82 model.
%
%    Examples:
%
%      gmst=ut2gmst('2012-01-04 15:00:03')
%      [gmst0,omegae]=ut2gmst('2012-01-04')
%      gmst=ut2gmst(datenum)

%   (c) Hans van der Marel, Delft University of Technology, 2012.
%
%   Created:    31 Dec 2012 by Hans van der Marel
%   Modified:   22 Jun 2013 by Hans van der Marel
%                 - minor changes to documentation

% Set default model

if nargin < 2 
    model='IAU-82';
end

% Convert input to Matlab date number

if ischar(ut1)
   ut=datenum(ut1); 
end

% GMST 

switch upper(model)
   case 'APPROXIMATE'
      t0=datenum([2000 1 1 12 0 0]);
      gmst = rem( 18.697374558 + 24.06570982441908 .* (ut1-t0) , 24);   % in hours
      gmst = gmst*3600;
   case 'IAU-82'
      % time since J2000 (Jan 1, 2000, 12 h) in Julian centuries
      j2000=(ut1 - datenum([2000 1 1 12 0 0])) ./ 36525.0;
      % gmst in seconds
      gmst = ( ( - 6.2e-6 .* j2000 + 0.093104 ) .* j2000   ...
         + ( 876600.0 * 3600.0 + 8640184.812866 ) ) .* j2000 + 67310.54841;
      gmst = mod(gmst, 86400); 
   otherwise
      error('ut2gmst unknown model')
end

% convert GMST to radians [0 2pi]

gmst = gmst.*pi./43200;

% Rotational velocity of the Earth [rev/day] 

% Me = 7.2921151467e-5;     [rad/sec]
omegae = 1.0027379093;

end


