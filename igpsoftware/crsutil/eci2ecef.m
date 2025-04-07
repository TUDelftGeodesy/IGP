function [xsate,vsate]=eci2ecef(t,xsat,vsat)
%ECI2ECEF    Convert position and velocity from ECI to ECEF reference frame.
%   [XSATE,VSATE]=ECI2ECEF(T,XSAT,VSAT) converts cartesian coordinates XSAT  
%   and velocities VSAT in Earth Centered Inertial (ECI) reference frame,
%   with T the time (UT1) given as Matlab datenumbers, into cartesian 
%   coordinates XSATE and velocities VSATE in Earth Centered Earth Fixed (ECEF) 
%   reference frame. T is a vector with length n, the number of epochs, and
%   XSAT and VSAT are n-by-3 matrices.
%
%   XSATE=ECI2ECEF(T,XSAT) only transforms the positions.  
%
%   See also UT2GMST and ECEF2ECI.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016.

%   Created:  26 November 2016 by Hans van der Marel
%   Modified: 

% Check input arguments

if nargin < 2 || nargin > 3
    error('Incorrect number of arguments.')
end
if nargin < 3
    vsat=nan(size(xsat));
end

if size(t,1) ~= size(xsat,1)
    error('Size of T does not match XSAT.')
end
if size(xsat,2)~=3
    error('XSAT incorrect size.')
end
if size(vsat,1) ~= size(xsat,1)
    error('Size of VSAT does not match XSAT.')
end
if size(vsat,2)~=3
    error('VSAT incorrect size.')
end

% Compute rotation angle (GMST) around Z-axis

[gst0,omegae]=ut2gmst(t(1));
gst=gst0+2*pi*omegae*(t-t(1));

% Rotate satellite positions around z-axis (ECI->ECEF))

xsate(:,1)=  cos(gst).*xsat(:,1) + sin(gst).*xsat(:,2);
xsate(:,2)= -sin(gst).*xsat(:,1) + cos(gst).*xsat(:,2);
xsate(:,3)=  xsat(:,3);

% To convert the velocity is more complicated. The velocity in ECEF
% consists of two parts. We find this by differentiating the transformation
% formula for the positions
%
%    xsate = R * xsat
%
% This gives (product rule, and some rewriting), with |_dot| the derivatives
%
%    xsate_dot = R * xsat_dot + R_dot * xsat    <=>
%    vsate = R * ( vsat + inv(R)*R_dot * xsat ) <=>
%    vsate = R * ( vsat + W * xsat )
%
% with |W = inv( R )*R_dot = [ 0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ]| and
% with |w| the angular velocity vector of the ECEF frame with respect to
% the ECU frame, expressed in the ECI frame.
% 
% The velocity vector in the ECEF is computed as follows

w=[0;0;2*pi*omegae/86400];
W = [ 0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0];

vsatw=vsat-xsat*W';

vsate(:,1)=  cos(gst).*vsatw(:,1) + sin(gst).*vsatw(:,2);
vsate(:,2)= -sin(gst).*vsatw(:,1) + cos(gst).*vsatw(:,2);
vsate(:,3)=  vsatw(:,3);

end