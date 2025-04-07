function [orb,corbtype]=vec2orb(svec,GM)
%VEC2ORB   Convert inertial state vector into Keplerian elements.
%   ORB=VEC2ORB(SVEC) converts a 6-element inertial state vector SVEC
%   with cartesian position and velocity (X, Y, Z, Xdot, Ydot, Zdot) into
%   a vector ORB with 6 Keplerian elements, with
%         ORB(:,1)    Semi-major axis (meters),
%         ORB(:,2)    Eccentricity (unity),
%         ORB(:,3)    Inclination (radians),
%         ORB(:,4)    Right ascension of the ascending node (radians),
%         ORB(:,5)    Argument of the pericenter (radians),
%         ORB(:,6)    True anomaly (radians).
%   This routine is fully vectorized, if SVEC is a matrix then ORB will be
%   a matrix of the same size. One of the dimensions of the input matrix must 
%   be 6. Units are meter, meter/sec or radians.  
%
%   ORB=VEC2ORB(SVEC,GM) provides an optional gravitational parameter of the 
%   central body. Default for GM [meter**3/sec**2] is the IERS 1996 standard  
%   value for the Earth (GM=3986004418e5)
%
%   [ORB,CORBTYPE]=VEC2ORB(SVEC) also outputs a character array CORBTYPE
%   with the orbit type, with
%       'ei'   elliptical inclined     (all Kepler elements defined)
%       'ci'   circular inclined       (w =0, nu=arglat)
%       'ee'   elliptical equatorial   (w=lonper, omega=0)
%       'ce'   circular equatorial     (w=0, omega=0, nu=truelon)
%   orbits are "circular" if the eccentricity < eps and "equatorial" if
%   the inclination < eps, with eps=1e-8.  
%
%   See also ORB2VEC, KEPLER, KEPLERINV

%   (c) Hans van der Marel, Delft University of Technology, 2010.
%
%   Created:    16 Nov 2010 by Hans van der Marel
%   Modified:   

% Check the input arguments

if nargin <2, GM=3986004418e5;, end      % Set default value for GM

[m,n]=size(svec);                        % Optionally transpose svec and check size
if ( m == 6 && n ~= 6 )
   transpose=true;
   svec=svec';
   [m,n]=size(svec);
else
   transpose=false;
end
if ( n ~= 6 ) 
   error('Input must be a vector/matrix with 6 state vector elements')
end

% Inner products (rrdot = R.V , r=sqrt(R.R) , vsq = V.V )

rrdot=svec(:,1).*svec(:,4)+svec(:,2).*svec(:,5)+svec(:,3).*svec(:,6);
r=sqrt(svec(:,1).*svec(:,1)+svec(:,2).*svec(:,2)+svec(:,3).*svec(:,3));
vsq=svec(:,4).*svec(:,4)+svec(:,5).*svec(:,5)+svec(:,6).*svec(:,6);

% Angular momentum vector (H = R x V)

hx=svec(:,2).*svec(:,6)-svec(:,3).*svec(:,5);
hy=svec(:,3).*svec(:,4)-svec(:,1).*svec(:,6);
hz=svec(:,1).*svec(:,5)-svec(:,2).*svec(:,4);

hsini2=hx.*hx+hy.*hy;
hsq=hsini2+hz.*hz;
h=sqrt(hsq);

% Semi-major axis

ainv=2./r-vsq./GM;
a=1./ainv;

% Eccentricity

ome2=hsq.*ainv./GM;
ecc=sqrt(1.0 - ome2);
ecc(ome2 > 1)=0;                  % special handling of negative values

% Inclination (0...pi)

incl=acos(hz./h);

% Determine orbit type (for handling of special cases)

small=1e-8;

idxecc=(ecc < small);
idxincl=( (incl < small) | (abs(incl-pi) < small ));

idx_ce = ( idxecc & idxincl );      % circular equatorial => w=0, omega=0, nu=truelon
idx_ci = ( idxecc & ~idxincl );     % circular inclined => w =0, nu=arglat
idx_ee = ( ~idxecc & idxincl );     % elliptical equatorial => w=lonper, omega=0
idx_ei = ( ~idxecc & ~idxincl );    % elliptical inclined

orbtype(idx_ei)=0;
orbtype(idx_ee)=1;
orbtype(idx_ci)=2;
orbtype(idx_ce)=3;

corbdef=['ei';'ee';'ci';'ce'];
corbtype=corbdef(orbtype+1,:);

% Standard handling of elliptical inclined orbits...

% The computations below do not do special hanling of circular or equatorial
% orbits. This is possible because atan2(0,0)=0 is defined in Matlab, however,
% the some of the angles will actually be undefined.

% Longitude of ascending node (0...2*pi)

omega=atan2(hx,-hy);    
idx=find(omega < 0);
omega(idx)=omega(idx)+2*pi;

% True anomaly (0...2*pi)

resinf=a.*ome2.*rrdot./h;
recosf=a.*ome2-r;
nu=atan2(resinf,recosf);
idx=find(nu < 0);
nu(idx)=nu(idx)+2*pi;

% Argument of perigee (0...2*pi)

suprod=-hz.*(svec(:,1).*hx+svec(:,2).*hy)+svec(:,3).*hsini2;
cuprod=h.*(-svec(:,1).*hy+svec(:,2).*hx);
w=atan2(suprod.*recosf-cuprod.*resinf,cuprod.*recosf+suprod.*resinf);
idx=find(w < 0);
w(idx)=w(idx)+2*pi;

% Special handling of circular and equatorial orbits

% not yet implemented

% Done, store Keplerian elements in ORB

orb=[ a ecc incl omega w nu];

if transpose                            
   orb=orb';
end

return
