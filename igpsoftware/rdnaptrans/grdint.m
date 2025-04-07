function zint=grdint(gstruct,xint,yint)
%GRDINT  Grid interpolation using Catmull-Rom-Overhauser splines
%  Grid interpolation using Catmull-Rom-Overhauser splines:
%  Input:
%    gstruct    structure with grid data
%      gstruct.x    vector with the x coordinates (longitude)
%      gstruct.y    vector with the y coordinates (latitude)
%      gstruct.mat  matrix with values of dimension [length(x),length(y)]
%    xint       x coordinate (longitude) of the interpolated point
%    yint       y coordinate (lattitude) of the interpolated point
%  Output:
%    zint       interpolated values
%
%  Layout of mat:
%
%    +--> y                                        West 
%    |                                              ^
%    v   (1,1)  (1,2) .. (1,ny)                     |
%        (2,1)  (2,2)                       South <-+-> North
%    x     :                                        |
%        (nx,1)          (nx,ny)                    v
%                                                  East
%
%  See also readgrd and testgrdint.
%
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013.

%  Created:   6 Dec 2004 by Hans van der Marel, TUD
%  Modified  23 May 2011 by Piers Titus van der Torren
%              - Vectorized the function (initial version)
%             6 Jun 2013 by Hans van der Marel
%              - Added some corrections to the vectorization
%              - Updated comments

x=gstruct.x;
y=gstruct.y;
mat=gstruct.mat;

% Check the size of the vectors and arrays

[nx,ny]=size(mat);
if nx ~= length(x)
   error('length of array with x values does not match matrix');
end
if ny ~= length(y)
   error('length of array with y values does not match matrix');
end

% Check if the point to be interpolated is within the grid bounding box

%invalid= xint < x(1) | xint > x(nx) | yint < y(1) | yint > y(ny) | isnan(xint) | isnan(yint);
%invalid= xint <= x(1) | xint >= x(nx) | yint <= y(1) | yint >= y(ny) | isnan(xint) | isnan(yint);
invalid= xint < x(2) | xint > x(nx-1) | yint < y(2) | yint > y(ny-1) | isnan(xint) | isnan(yint);

% Compute the step-size

dx=(x(nx)-x(1))/(nx-1);
dy=(y(ny)-y(1))/(ny-1);


% The selected grid points are situated around point X like this:
%
%    +-->  y                                              
%    |                                                                  West
%    v   (ix-1,iy-1)  (ix-1,iy)    (ix-1,iy+1)   (ix-1,iy+2)             ^
%                                                                        |
%    x   (ix,iy-1)    (ix,iy)      (ix,iy+1)     (ix,iy+2)       South <-+-> North
%                               X                                        |
%        (ix+1,iy-1)  (ix+1,iy)    (ix+1,iy+1)   (ix+1,iy+1)             v
%                                                                       East
%        (ix+2,iy-1)  (ix+2,iy)    (ix+2,iy+1)   (ix+2,iy+2)             
%                                                         
%        
% Find the grid value just to the left and above the target point

ix=floor((xint-x(1))/dx)+1;
iy=floor((yint-y(1))/dy)+1;

ix(invalid)=0;
iy(invalid)=0;

ix=max(ix,2);
ix=min(ix,nx-2);
iy=max(iy,2);
iy=min(iy,ny-2);

ddx = (xint-x(1))/dx-ix+1;
ddy = (yint-y(1))/dy-iy+1;

% Compute weighting factors for x and y components

ddx2=ddx.^2;
ddx3=ddx.^3;
ddy2=ddy.^2;
ddy3=ddy.^3;
f= [   -0.5*ddx+ddx2-0.5*ddx3, ...
        1.0-2.5*ddx2+1.5*ddx3, ...
    0.5*ddx+2.0*ddx2-1.5*ddx3, ...
           -0.5*ddx2+0.5*ddx3];
g= [   -0.5*ddy+ddy2-0.5*ddy3, ...
        1.0-2.5*ddy2+1.5*ddy3, ...
    0.5*ddy+2.0*ddy2-1.5*ddy3, ...
           -0.5*ddy2+0.5*ddy3];
    
% Compute interpolated value ( zint=f'*mat(ix-1:ix+2,iy-1:iy+2)*g )

idx = sub2ind([nx,ny],ix,iy);
zint=zeros(size(xint));
for lx=-1:2
   for ly=-1:2
      zint = zint + f(:,lx+2) .* g(:,ly+2) .* mat(idx +lx + ly*nx);
   end
end

% Replace out of bound points with NaN

zint(invalid)=nan;

return;
