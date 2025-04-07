function [M,E]=keplernu(nu,ecc)
%KEPLERNU  Compute mean anomaly from true anomaly (Kepler's equation).
%   M=KEPLERNU(NU,ECC) computes the mean anomaly M [rad] from the true 
%   anomaly NU [rad] and eccentricity ECC. This routine is fully vectorized, 
%   if NU is a vector, M will be a vector with the same dimensions.
%
%   [M,E]=KEPLERNU(NU,ECC) returns also the eccentric anomaly E [rad].
%
%   This routine should only be used for elliptical orbits. Parabolic and 
%   hyperbolic orbits are not supported and give false results (this is
%   nowhere checked for).
%
%   See also KEPLERM and KEPLER.

%   (c) Hans van der Marel, Delft University of Technology, 2012.
%
%   Created:    31 Dec 2012 by Hans van der Marel
%   Modified:   22 Jun 2013 by Hans van der Marel
%                 - minor changes to documentation

% Compute eccentric anomaly

denom = 1.0 + ecc.*cos(nu); 

sine= ( sqrt( 1.0 -ecc.*ecc ) .* sin(nu) ) ./ denom ;
cose= ( ecc + cos(nu) ) ./ denom;
E  = atan2( sine,cose );

% Compute mean anomaly

M = E - ecc.*sin(E);

return
