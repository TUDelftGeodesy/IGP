function [xsat,vsat]=ecef2eci(t,xsate,vsate)
%ECEF2ECI    Convert position and velocity from ECEF to ECI reference frame.
%   [XSAT,VSAT,]=ECEF2ECI(T,XSATE,VSATE) converts cartesian coordinates XSATE  
%   and velocities VSATE in Earth Centered Earth Fixed (ECEF) reference frame,
%   with T the time (UT1) given as Matlab datenumbers, into cartesian 
%   coordinates XSAT and velocities VSAT in Earth Centered Inertial (ECI) 
%   reference frame. T is a vector with length n, the number of epochs, and
%   XSATE and VSATE are n-by-3 matrices.
%
%   [XSAT,VSAT]=ECEF2ECI(T,XSATE) assumes the velocity in ECEF is zero.
%
%   See also UT2GMST and ECI2ECEF.
%
%   (c) Hans van der Marel, Delft University of Technology, 2016.

%   Created:  26 November 2016 by Hans van der Marel
%   Modified: 

% Check input arguments

if nargin < 2 || nargin > 3
    error('Incorrect number of arguments.')
end
if size(xsate(:),1) == 3 && size(t,1) > 1
    xsate=repmat(xsate(:)',size(t));
end

if size(t,1) ~= size(xsate,1)
    error('Size of T does not match XSATE.')
end
if size(xsate,2)~=3
    error('XSATE incorrect size.')
end
if size(t,1) ~= size(xsate,1)
    error('Size of T does not match XSATE.')
end

if nargin < 3
    vsate=zeros(size(xsate));
end

if size(vsate,1) ~= size(xsate,1)
    error('Size of VSATE does not match XSATE.')
end
if size(vsate,2)~=3
    error('VSATE incorrect size.')
end

% Compute rotation angle (GMST) around Z-axis

[gst0,omegae]=ut2gmst(t(1));
gst=gst0+2*pi*omegae*(t-t(1));
gst=-gst;

% Rotate satellite positions around z-axis (ECEF->ECI))

xsat(:,1)=  cos(gst).*xsate(:,1) + sin(gst).*xsate(:,2);
xsat(:,2)= -sin(gst).*xsate(:,1) + cos(gst).*xsate(:,2);
xsat(:,3)=  xsate(:,3);

% To convert the velocity is more complicated. The velocity in ECEF
% consists of two parts. We find this by differentiating the transformation
% formula for the positions
%
%    xsat = R * xsate
%
% This gives (product rule, and some rewriting), with |_dot| the derivatives
%
%    xsat_dot = R * xsate_dot + R_dot * xsate    <=>
%    vsat = R * ( vsate + inv(R)*R_dot * xsate ) <=>
%    vsat = R * ( vsate + W * xsate )
%
% with |W = inv( R )*R_dot = [ 0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0 ]| and
% with |w| the angular velocity vector of the ECI frame with respect to
% the ECEF frame, expressed in the ECEF frame.
% 
% The velocity vector in the ECI is computed as follows

w=[0;0;-2*pi*omegae/86400];
W = [ 0 -w(3) w(2) ; w(3) 0 -w(1) ; -w(2) w(1) 0];

vsatw=vsate-xsate*W';

vsat(:,1)=  cos(gst).*vsatw(:,1) + sin(gst).*vsatw(:,2);
vsat(:,2)= -sin(gst).*vsatw(:,1) + cos(gst).*vsatw(:,2);
vsat(:,3)=  vsatw(:,3);

end