function testgrdint(grdfile)
%TESTGRDINT  Test grid interpolation function for RDNAPTRANS
%  This funnction tests the grid interpolation for RDNAPTRANS.
%  Syntax
%             testgrdint()
%             testgrdint('nlgeo04.grd')
%             testgrdint('x2c.grd')
%             testgrdint('y2c.grd')
%
%  See also readgrd and grdint.
% 
%  (c) Hans van der Marel, Delft University of Technology, 2004-2013

% Created:   6 Dec 2004 by Hans van der Marel, TUD
% Modified: 11 Jun 2013 by Hans van der Marel, TUD
%             - converted into function
%             - test vectorized version of grid interpolation

% Matlab script to test the grid interpolation

if nargin < 1
   grdfile='nlgeo04.grd';
end

% If scalar =  1 test also scalar calls (will be time consuming)

scalar=0;

% Read the grid 

NLGRD=readgrd(grdfile);

x=NLGRD.x;
y=NLGRD.y;

figure;imagesc(NLGRD.mat');colorbar;set(gca,'Ydir','normal');title('original grid');

% Test the interpolation on the grid values

[yi,xi]=meshgrid(y,x);
test3=grdint(NLGRD,xi(:),yi(:));
test3=reshape(test3,size(NLGRD.mat));
maxdiff=max(abs(NLGRD.mat(:)-test3(:)));
fprintf('max diff at grid values is       %9.5f\n', maxdiff)
if maxdiff > 1e-6
   figure;imagesc(NLGRD.mat'-test3');colorbar;set(gca,'Ydir','normal');title('difference with interpolation');
end

% Increase the grid samples by a factor 5

ix=1;
iy=1;
nx=(length(x)-1)*5+1;
ny=(length(y)-1)*5+1;
dx=(x(ix+1)-x(ix))/5;
dy=(y(iy+1)-y(iy))/5;

[yi,xi]=meshgrid(y(iy)+[0:ny-1]'*dy,x(ix)+[0:nx-1]'*dx);

% test GRDINT in vector mode mode

tstart=tic;
test2=grdint(NLGRD,xi(:),yi(:));
test2=reshape(test2,[nx ny]);
fprintf('vectorized call, elapsed time is %9.2f\n', toc(tstart) );
figure;imagesc(test2');colorbar;set(gca,'Ydir','normal');title('interpolated grid (vectorized call)');

% test GRDINT in scalar mode

if scalar
  test1=zeros(ny*nx,1);
  tstart=tic;
  for l=1:ny*nx
    test1(l)=grdint(NLGRD,xi(l),yi(l));      
  end
  test1=reshape(test1,[nx ny]);
  fprintf('scalar call, elapsed time is     %9.2f\n', toc(tstart) );
  figure;imagesc(test1');colorbar;set(gca,'Ydir','normal');title('interpolated grid (scalar calls)');
  fprintf('max diff scalar vs vectorized is %9.5f\n', max(abs(test1(:)-test2(:))))
end

return;
