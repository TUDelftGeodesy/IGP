function str=plh2str(plh,fmt)
%PLH2STR Convert ellipsoidal coordinates phi,lambda,h to string.
%   STR=PLH2STR(PLH) converts the Nx3 matrix PLH with in the rows ellipsoidal
%   coordinates Phi, Lambda and Height into a character string with format
%   FMT='  %3d %02d %07.4f %4d %02d %07.4f  %9.4f' (default). Phi and Lambda 
%   are in radians, h is in meters.
%
%   STR=PLH2STR(PLH,FMT) uses the format FMT instead of the default format.
%   If FMT='pretty' it prints degree, minute and second symbols in the
%   default format.
%
%   See also STR2PLH, PLH2XYZ and XYZ2PLH.
%
%   (c) Hans van der Marel, Delft University of Technology, 1995,2013

%   Created:    29 Apr 1995 by Hans van der Marel
%   Modified:   14 Jun 2013 by Hans van der Marel
%                - updated description
%               12 March 2014 by Hans van der Marel
%                - fix for Southern and Western hemispheres
%                - pretty print option

pretty=false;
if nargin==1
    fmt='  %3d %02d %07.4f %4d %02d %07.4f  %9.4f';
elseif nargin ==2
    if strcmpi(fmt,'pretty')
       pretty=true;
       fmt='  %3d%c%02d''%07.4f" %4d%c%02d''%07.4f"  %9.4f';
    end
else
    error('incorrect number of input arguments')
end

rad=45/atan(1);
phi(:,1)=abs(plh(:,1))*rad;
lam(:,1)=abs(plh(:,2))*rad;

phi(:,3)=rem(phi(:,1)*3600,60);
phi(:,2)=fix(rem(phi(:,1)*60,60));
phi(:,1)=sign(plh(:,1)).*fix(phi(:,1));

lam(:,3)=rem(lam(:,1)*3600,60);
lam(:,2)=fix(rem(lam(:,1)*60,60));
lam(:,1)=sign(plh(:,2)).*fix(lam(:,1));

for i=1:size(plh,1)
  if pretty
     str(i,:)=sprintf(fmt,phi(i,1),char(176),phi(i,2),phi(i,3),lam(i,1),char(176),lam(i,2),lam(i,3),plh(i,3));
  else
     str(i,:)=sprintf(fmt,phi(i,1),phi(i,2),phi(i,3),lam(i,1),lam(i,2),lam(i,3),plh(i,3));
  end
end

return