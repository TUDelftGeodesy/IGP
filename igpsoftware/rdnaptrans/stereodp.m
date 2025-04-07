function [out1,out2]=stereodp(in1,in2,mstruct,direction)
%STEREODP Stereographic double projection 
%  Stereographic double projection
%  Syntax
%                [x,y]=stereodp(phi,lambda,mstruct,'forward')
%                [phi,lambda]=stereodp(x,y,mstruct,'inverse')
%  Input
%    mstruct     structure with map-projection data
%    direction   direction of the projection, forward or reverse
%  Input/output
%    phi         latitude in radians
%    lambda      longitude in radians
%    x           x coordinate (easting)
%    y           y coordinate (norting)
%
%  References:
%    G. Bakker, J.C. de Munck and G.L. Strang van Hees, "Radio Positioning 
%      at Sea". Delft University of Technology, 1995.
%    G. Strang van Hees, "Globale en lokale geodetische systemen", 
%      Nederlandse Commissie voor Geodesie (NCG), Delft, 1997.
%
%  See also nlbessel.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013.

%  Created:   6 Dec 2004 by Hans van der Marel, TUD
%  Modified  23 May 2011 by Piers Titus van der Torren
%              - Vectorized the function
%             6 Jun 2013 by Hans van der Marel
%              - Updated comments and check for precision

PRECISION     = 0.0001;
RAD_PRECISION = PRECISION/40e6*pi*2;

% Extract the necessary projection data and convert to radians

phic  = mstruct.origin(1)*pi/180;
lamc  = mstruct.origin(2)*pi/180;
a     = mstruct.a;
inv_f = mstruct.inv_f;

scale = mstruct.scale;
falsex= mstruct.falseeasting;
falsey= mstruct.falsenorthing;

% Compute some more constants for use later on

  %   f              flattening of the ellipsoid
  %   ee             first eccentricity squared (e squared in some notations)
  %   e              first eccentricity 
  %   eea            second eccentricity squared (e' squared in some notations)

  %   phi0_sphere    latitude of projection origin on sphere
  %   lam0_sphere    longitude of projection origin on sphere

  %   r1             first (North South) principal radius of curvature in origin (M in some notations)
  %   r2             second (East West) principal radius of curvature in origin (N in some notations)
  %   r_sphere       radius of sphere

  %   n              constant of Gaussian projection n = 1.000475... 
  %   q0             isometric latitude of origin on ellipsiod
  %   w0             isometric latitude of origin on sphere
  %   m              constant of Gaussian projection m = 0.003773... (also named c in some notations)
 
f=1/inv_f;
ee=f*(2-f);  
e=sqrt(ee);  
eea = ee/(1.0-ee); 

phi0_sphere = atan(tan(phic)/sqrt(1+eea*cos(phic)^2));
lam0_sphere = lamc;

r1 = a*(1-ee)/sqrt(1-ee*sin(phic)^2)^3;
r2 = a/sqrt(1.0-ee*sin(phic)^2);
r_sphere = sqrt(r1*r2);

n = sqrt(1+eea*cos(phic)^4);
q0 = atanh(sin(phic))-e*atanh(e*sin(phic));
w0 = log(tan(pi/4+0.5*phi0_sphere));
m = w0-n*q0; 

% Do the actual projection and inverse projection

switch direction
case 'forward'

  % Forward projection

  phi=in1;
  lambda=in2;

  % q                    isometric latitude on ellipsiod
  % w                    isometric latitude on sphere
  % phi_sphere           latitide on sphere
  % delta_lambda_sphere  difference in longitude on sphere with origin
  % psi                  distance angle from origin on sphere
  % alpha                azimuth from origin
  % r                    distance from origin in projection plane

  q                    = atanh(sin(phi))-e*atanh(e*sin(phi));
  w                    = n*q+m;
  phi_sphere           = 2*atan(exp(w))-pi/2;
  delta_lambda_sphere  = n*(lambda-lam0_sphere);
  sin_half_psi_squared = sin(0.5*(phi_sphere-phi0_sphere)).^2 + sin(0.5*delta_lambda_sphere).^2 .* cos(phi_sphere).*cos(phi0_sphere);
  sin_half_psi         = sqrt(sin_half_psi_squared);
  cos_half_psi         = sqrt(1-sin_half_psi_squared);
  tan_half_psi         = sin_half_psi./cos_half_psi;
  sin_psi              = 2*sin_half_psi.*cos_half_psi;
  cos_psi              = 1-2*sin_half_psi_squared;
  sin_alpha            = sin(delta_lambda_sphere).*(cos(phi_sphere)./sin_psi);
  cos_alpha            = (sin(phi_sphere)-sin(phi0_sphere).*cos_psi)./(cos(phi0_sphere).*sin_psi);
  r                    = 2*scale*r_sphere.*tan_half_psi;

  %  Set the output variables

  out1 = r.*sin_alpha+falsex;
  out2 = r.*cos_alpha+falsey; 

case 'inverse'

  % Inverse projection

  x=in1;
  y=in2;

  % r                    distance from origin in projection plane
  % alpha                azimuth from origin
  % psi                  distance angle from origin on sphere
  % phi_sphere           latitide on sphere 
  % delta_lambda_sphere  difference in longitude on sphere with origin
  % w                    isometric latitude on sphere
  % q                    isometric latitude on ellipsiod

  r  = sqrt((x-falsex).^2+(y-falsey).^2);
  issmall=(r < PRECISION);
  rr=r;
  rr(issmall)=PRECISION;
  sin_alpha = (x-falsex)./rr;
  cos_alpha = (y-falsey)./rr;
  sin_alpha(issmall) = 0;
  cos_alpha(issmall) = 1;
  psi                 = 2*atan(r/(2*scale*r_sphere));
  phi_sphere          = asin(cos_alpha.*cos(phi0_sphere).*sin(psi)+sin(phi0_sphere).*cos(psi));
  delta_lambda_sphere = asin((sin_alpha.*sin(psi))./cos(phi_sphere));

  lambda = delta_lambda_sphere/n+lamc;

  w = atanh(sin(phi_sphere));
  q = (w-m)/n;
    
  %  Iterative calculation of phi
  
  phi=0;
  diff=pi/2;
  while any(diff > RAD_PRECISION)
    previous = phi;
    phi      = 2*atan(exp(q+0.5*e*log((1+e*sin(phi))./(1-e*sin(phi)))))-pi/2;
    diff     = abs(phi-previous);
  end

  %  Set the output variables

  out1 = phi;
  out2 = lambda; 

otherwise

  error('Unrecognized direction string')

end

return;


