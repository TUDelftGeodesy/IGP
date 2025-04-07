function M=kepler(E,ecc)
%KEPLER    Compute mean anomaly from eccentric anomaly (Kepler's equation).
%   M=KEPLER(E,ECC) computes the mean anomaly M [rad] from the eccentric 
%   anomaly E [rad] and eccentricity ECC. This routine is fully vectorized, 
%   if E is a vector, M will be a vector with the same dimensions.
%
%   This routine should only be used for elliptical orbits. Parabolic and 
%   hyperbolic orbits are not supported and give false results (this is
%   nowhere checked for).
%
%   See also KEPLERM and KEPLERNU.

%   (c) Hans van der Marel, Delft University of Technology, 2010-2012.
%
%   Created:    20 Aug 2010 by Hans van der Marel
%   Modified:   22 Jun 2013 by Hans van der Marel
%                 - minor changes to documentation

M = E - ecc.*sin(E);

return
