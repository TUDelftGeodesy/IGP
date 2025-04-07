function plh=str2plh(str,fmt)
%STR2PLH Read ellipsoidal coordinates phi,lambda,h from string.
%   PLH=STR2PLH(STR) converts the character array STR with ellipsoidal coordi-
%   nates into a matrix PLH with in the rows ellipsoidal coordinates Phi, 
%   Lambda and Height. Phi and Lambda are converted to radians, h is in meters.
%
%   PLH=STR2PLH(STR,FMT) uses the format FMT instead of the default format
%   FMT='%f %f %f %f %f %f  %f'.
%
%   See also STR2PLH, PLH2XYZ and XYZ2PLH.
%
%   (c) Hans van der Marel, Delft University of Technology, 2001,2013

%   Created:    30 Oct 2001 by Hans van der Marel
%   Modified:   14 Jun 2013 by Hans van der Marel
%                - updated description
%               12 March 2014 by Hans van der Marel
%                - fix for Southern and Western hemispheres
%                - strip characters from pretty print option

if nargin==1, fmt='%f %f %f %f %f %f  %f';, end 

rad=45/atan(1);

tmp=zeros(size(str,1),7);
for i=1:size(str,1)
   tmpstr=str(i,:);
   tmpstr=strrep(tmpstr,char(176),' ');
   tmpstr=strrep(tmpstr,'''',' ');
   tmpstr=strrep(tmpstr,'"',' ');
   tmp(i,:)=sscanf(tmpstr,fmt,7)';
end

plh(:,1)=sign(tmp(:,1)).*( abs(tmp(:,1))+tmp(:,2)./60+tmp(:,3)./3600 )/rad;
plh(:,2)=sign(tmp(:,4)).*( abs(tmp(:,4))+tmp(:,5)./60+tmp(:,6)./3600 )/rad;
plh(:,3)=tmp(:,7);

return