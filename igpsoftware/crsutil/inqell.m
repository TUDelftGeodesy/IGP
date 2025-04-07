function [a,f,GM] = inqell(ellips)
%INQELL  Semi-major axis, flattening and GM for various ellipsoids.
%   [A,F,GM]=INQELL(ELLIPS) returns the semi-major axis A, flattening F and GM 
%   for a particular ellipsoid ELLIPS. Default for ELLIPS is 'WGS-84'. If 
%   ELLIPS is '?' a list of supported ellipsoids is printed.
%
%   See also PLH2XYZ and XYZ2PLH.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995,2013

%   Created:     7 May 1995 by Hans van der Marel
%   Modified:   14 Jun 2013 by Hans van der Marel
%                - updated description

ell= ['AIRY         ';
      'BESSEL       ';
      'CLARKE       ';
      'INTERNATIONAL';
      'HAYFORD      ';
      'GRS80        ';
      'WGS-84       '];
par= [6377563.396 , 299.324964    , NaN        ;
      6377397.155 , 299.1528128   , NaN        ;
      6378249.145 , 293.465       , NaN        ;
      6378388.      297.00        , NaN        ;
      6378388.      297.00        , 3.986329e14;      
      6378137.    , 298.257222101 , 3.986005e14;
      6378137.    , 298.257223563 , 3.986005e14];

if nargin==0,ellips='unknown';, end
ellips=deblank(upper(ellips));
if strcmp(ellips,'?')
  disp('Ellipsoids:');
  disp(' ');
  disp(ell);
  return
end

i=0;
for j=1:size(par,1)
  if strcmp(deblank(ell(j,:)),ellips), i=j;, end
end
if i==0 
  i=size(par,1); 
  disp(['Warning: Ellipsoid ',ellips,' not found, ',deblank(ell(i,:)), ...
  ' selected instead']);
end

a=par(i,1);
f=1/par(i,2);
GM=par(i,3);
