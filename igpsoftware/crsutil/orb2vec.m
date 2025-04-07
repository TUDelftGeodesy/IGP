function svec=orb2vec(orb,GM)
%ORB2VEC   Convert Keplerian elements into inertial state vector.
%   SVEC=ORB2VEC(ORB) converts the vector ORB with 6 Keplerian elements
%   into the 6-element inertial state vector SVEC with cartesian position and 
%   velocity (X, Y, Z, Xdot, Ydot, Zdot). The Keplerian elements are
%         ORB(:,1)    Semi-major axis (meters),
%         ORB(:,2)    Eccentricity (unity),
%         ORB(:,3)    Inclination (radians),
%         ORB(:,4)    Right ascension of the ascending node (radians),
%         ORB(:,5)    Argument of the pericenter (radians),
%         ORB(:,6)    True anomaly (radians).
%   This routine is fully vectorized, if ORB is a matrix then SVEC 
%   will be a matrix of the same size. One of the dimensions of the input
%   matrix must be 6. The units are meter, meter/sec or radians.  
%
%   SVEC=ORB2VEC(ORB,GM) provides an optional gravitational parameter of the 
%   central body. Default for GM [meter**3/sec**2] is the IERS 1996 standard  
%   value for the Earth (GM=3986004418e5)
%
%   See also VEC2ORB, KEPLER, KEPLERINV

%   (c) Hans van der Marel, Delft university of Technology 2010.
%
%   Created:    16 Nov 2010 by Hans van der Marel
%   Modified:   

% Check the input arguments

if nargin <2, GM=3986004418e5;, end     % Set default value for GM

[m,n]=size(orb);                        % Optionally transpose orb and check size
if ( m == 6 && n ~= 6 )
   transpose=true;
   orb=orb';
   [m,n]=size(orb);
else
   transpose=false;
end
if ( n ~= 6 ) 
   error('Input must be a vector/matrix with 6 Keplerian elements')
end

% Compute position (rx,ry) and velocity (vx,vy) in orbital plane (perifocal system)

ecc = orb(:,2);                         % Eccentricity
cosnu = cos(orb(:,6));                  % Cosine and sine of true anomaly (nu)
sinnu = sin(orb(:,6));

p = orb(:,1).*( 1.0 - ecc.*ecc);        % Parameter of the ellipse p=a*(1-e^2) 

r = p ./ ( 1.0 + ecc.*cosnu);           % Length of position vector 

rx = r.*cosnu;                          % Position (rx,ry) in orbital plane
ry = r.*sinnu;

p(abs(p) < 0.0001)=0.0001;              % Protect against division by zero
tmp=sqrt(GM)./sqrt(p);

vx = -tmp.*sinnu;                       % Velocity (vx,vy) in orbital plane
vy = tmp.*(ecc + cosnu);

% Convert into inertial frame (3-1-3 Euler rotations)

cosincl=cos(orb(:,3));                  % Cosine and sine of inclination (incl)
sinincl=sin(orb(:,3));
cosomega=cos(orb(:,4));                 % Cosine and sine of longitude of ascending node (omega)
sinomega=sin(orb(:,4)); 
cosw=cos(orb(:,5));                     % Cosine and sine of argument of perigee (w)
sinw=sin(orb(:,5));

rx0 = cosw.*rx - sinw.*ry;              % Cosine and sine of argument of latitude u=w+nu          
ry0 = cosw.*ry + sinw.*rx;

vx0 = cosw.*vx - sinw.*vy;
vy0 = cosw.*vy + sinw.*vx;

svec= [ rx0.*cosomega  - ry0.*cosincl.*sinomega ...
        rx0.*sinomega  + ry0.*cosincl.*cosomega ...
                         ry0.*sinincl           ...
        vx0.*cosomega  - vy0.*cosincl.*sinomega ...
        vx0.*sinomega  + vy0.*cosincl.*cosomega ...
                         vy0.*sinincl           ];

% Let size of output vector/matrix match input

if transpose                            
   svec=svec';
end

return
