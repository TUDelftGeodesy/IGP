function [E,nu]=keplerm(M,ecc,TOL)
%KEPLERM  Compute eccentric/true from mean anomaly solving Kepler's eqn. 
%   E=KEPLERM(M,ECC) computes the eccentric anomaly E [rad] from the mean 
%   anomaly M [rad] and eccentricity ECC by solving Kepler's equation 
%   M=E-ECC*sin(E) iteratively using Newton's method. This routine is fully 
%   vectorized, if M is a vector, E will be a vector with the same dimensions.
%
%   [E,NU]=KEPLERM(M,ECC) returns also the true anomaly NU [rad].
%
%   E=KEPLERM(M,ECC,TOL) uses TOL as the cutoff criterion for the iterations, 
%   the default value is 1e-10.
%
%   This routine should only be used for elliptical orbits. Parabolic and 
%   hyperbolic orbits are not supported and give false results (this is
%   nowhere checked for).
%
%   See also KEPLERNU and KEPLER.

%   (c) Hans van der Marel, Delft University of Technology, 2010-2012.
%
%   Created:    20 Aug 2010 by Hans van der Marel
%   Modified:   31 Dec 2012 by Hans van der Marel
%               - renamed the routine from keplerinv to keplerm
%               - added computation of true anomaly
%               22 Jun 2013 by Hans van der Marel
%                 - minor changes to documentation

if nargin < 3, TOL=1e-10; end    % Set value for the tollerance

E = M;                           % Use M for the first value of E

f = ones(size(E));               % Newton's method for root finding
while max(abs(f)) > TOL
   f = M - E + ecc.*sin(E);      % Kepler's Equation
   fdot = -1 + ecc.*cos(E);      % Derivative of Kepler's equation
   E = E - f./fdot;
end

if nargout < 2, return, end
    
% sinv= ( sqrt( 1.0 -ecc*ecc ) * sin(E) ) / ( 1.0 -ecc*cos(E) );
% cosv= ( cos(E)-ecc ) / ( 1.0  - ecc*cos(E) );

sinnu = -1 * sqrt( 1 - ecc.*ecc ) .* sin(E) ./ fdot;
cosnu = ( ecc - cos(E) ) ./ fdot;

nu = atan2( sinnu,cosnu );       % True anomaly

return
