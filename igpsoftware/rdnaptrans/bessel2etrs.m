function trf=bessel2etrs(id)
%BESSEL2ETRS Transformation parameters for NL-Bessel to ETRS89 datum  
%  Transformation parameters for NL-Bessel to ETRS89 datum  .
%  Syntax
%          trf=bessel2etrs
%          trf=bessel2etrs(id)
%  Input
%    id      char       system id RDNAPTRANS2004 or RDNAPTRANS2008 (default 
%                       most recent
%  Output
%    trf   Structure with transformation parameters
%            trf.tx     translation in direction of x axis in meters
%            trf.ty     translation in direction of y axis in meters
%            trf.tz     translation in direction of z axis in meters
%            trf.alpha  rotation around x axis in radians
%            trf.beta   rotation around y axis in radians
%            trf.gamma  rotation around z axis in radians
%            trf.delta  scale parameter (scale = 1 + delta)
%            trf.x0     x-coordinate of pivot point
%            trf.y0     y-coordinate of pivot point
%            trf.z0     z-coordinate of pivot pount
%
%  Transformation parameters are relative to a pivot point ! The pivot
%  point is Amersfoort.
%
%  See also simtrans and etrs2bessel.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013

%  Created:   6 Dec 2004 by Hans van der Marel, TUD
%  Modified: 22 Apr 2010 by Sierd de Vries
%               - Updated parameters for RDNAPTRANS2008
%             6 Jun 2013 by Hans van der Marel
%               - Select between RDNAPTRANS2004 and RDNAPTRANS2008 

if (nargin < 1) 
   id='MOSTRECENT';
end

switch upper(id)
   case {'RDNAPTRANS2008','2008','MOSTRECENT','CURRENT'}

     trf.id    =  'RDNAPTRANS2008';

     % Transformation parameters from RD(Bessel) to ETRS89
     % Pivot point: center of the ellipsoid
     % tx           565.4171      m
     % ty            50.3319      m
     % tz           465.5524      m
     % alfa           1.9342      *10^-6 rad
     % beta          -1.6677      *10^-6 rad
     % gamma          9.1019      *10^-6 rad
     % delta          4.0725      *10^-6

     % Pivot point: Amersfoort
     % tx           593.0248      m
     % ty            25.9984      m
     % tz           478.7459      m
     % alfa           1.9342      *10^-6 rad
     % beta          -1.6677      *10^-6 rad
     % gamma          9.1019      *10^-6 rad
     % delta          4.0725      *10^-6
     % Cartesian coordinates of Amersfoort in RD(Bessel)
     % x0  3903453.1482           m
     % y0  368135.3134            m
     % z0  5012970.3051           m

     trf.tx    =     593.0248;
     trf.ty    =      25.9984;
     trf.tz    =     478.7459;
     trf.alpha =    1.9342e-6;
     trf.beta  =   -1.6677e-6;
     trf.gamma =    9.1019e-6;
     trf.delta =    4.0725e-6;
     trf.x0    = 3903453.1482;
     trf.y0    =  368135.3134;
     trf.z0    = 5012970.3051;

   case {'RDNAPTRANS2004','2004'}

     trf.id    =  'RDNAPTRANS2004';

     % Transformation parameters from RD(Bessel) to ETRS89
     % Pivot point: center of the ellipsoid
     % tx       565.2369            m
     % ty        50.0087            m
     % tz       465.6580            m
     % alpha      1.9725            *10^-6 rad
     % beta      -1.7004            *10^-6 rad
     % gamma      9.0677            *10^-6 rad
     % delta      4.0812            *10^-6

     % Transformation parameters from RD(Bessel) to ETRS89
     % Pivot point: Amersfoort
     % tx       593.0297            m
     % ty        26.0038            m
     % tz       478.7534            m
     % alpha      1.9725            *10^-6 rad
     % beta      -1.7004            *10^-6 rad
     % gamma      9.0677            *10^-6 rad
     % delta      4.0812            *10^-6
     % Cartesian coordinates of Amersfoort in RD(Bessel)
     % x0  3903453.1482           m
     % y0  368135.3134            m
     % z0  5012970.3051           m

     trf.tx    =     593.0297;
     trf.ty    =      26.0038;
     trf.tz    =     478.7534;
     trf.alpha =    1.9725e-6;
     trf.beta  =   -1.7004e-6;
     trf.gamma =    9.0677e-6;
     trf.delta =    4.0812e-6;
     trf.x0    = 3903453.1482;
     trf.y0    =  368135.3134;
     trf.z0    = 5012970.3051;

   otherwise
     error(['Unknown system ' id ])
end

return;
