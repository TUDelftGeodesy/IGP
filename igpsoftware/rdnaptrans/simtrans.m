function xyzout=simtrans(xyzin,trf)
%SIMTRANS  3D-simularity transformation for cartesian coordinates.
%  3D-simularity transformation for cartesian coordinates around a
%  pivot point using translation (3), rotation (3) and scale (1).
%  Syntax
%           xyzout=simtrans(xyzin,trf)
%  Input
%    xyzin  n-by-3 array with input coordinates [X,Y,Z]
%    trf    Structure with transformation parameters
%             trf.tx     translation in direction of x axis in meters
%             trf.ty     translation in direction of y axis in meters
%             trf.tz     translation in direction of z axis in meters
%             trf.alpha  rotation around x axis in radials
%             trf.beta   rotation around y axis in radials
%             trf.gamma  rotation around z axis in radials
%             trf.delta  scale parameter (scale = 1 + delta)
%             trf.x0     x-coordinate of pivot point
%             trf.y0     y-coordinate of pivot point
%             trf.z0     z-coordinate of pivot pount
%  Output
%    xyzout n-by-3 array with output coordinates [X,Y,Z]
%
%  See also bessel2etrs and etrs2bessel.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013.

%  Created:   6 Dec 2004 by Hans van der Marel, TUD
%  Modified  23 May 2011 by Piers Titus van der Torren
%              - Vectorized the function (initial version)
%             6 Jun 2013 by Hans van der Marel
%              - Updated vectorization part and harmonized input/output
%                arrays to n-by-3 types of arrays

% input must be n-by-3 array 

if all(size(xyzin) == [3,1])
  xyzin=xyzin';
  is3x1=true;
else
  is3x1=false;
end
assert(size(xyzin,2)==3,'input array must have shape n-by-3')

%  Compute the rotation matrix

alpha=trf.alpha;beta=trf.beta;gamma=trf.gamma;

R= [ cos(gamma)*cos(beta)                                    ... 
     cos(gamma)*sin(beta)*sin(alpha)+sin(gamma)*cos(alpha)   ...
    -cos(gamma)*sin(beta)*cos(alpha)+sin(gamma)*sin(alpha) ; ...
    -sin(gamma)*cos(beta)                                    ...
    -sin(gamma)*sin(beta)*sin(alpha)+cos(gamma)*cos(alpha)   ...
     sin(gamma)*sin(beta)*cos(alpha)+cos(gamma)*sin(alpha) ; ...
     sin(beta)                                               ...
    -cos(beta)*sin(alpha)                                    ...
     cos(beta)*cos(alpha)                                  ];

%  Apply the transformation

delta=trf.delta;
xyzpivot=[trf.x0,trf.y0,trf.z0];
xyztrans=[trf.tx,trf.ty,trf.tz];

xyzout = bsxfun(@plus, bsxfun(@minus,xyzin,xyzpivot)*R'.*(1.0+delta), xyztrans+xyzpivot);

if is3x1
   xyzout=xyzout';
end

return;


