function [lat,lon]=proj4inv(proj,x,y)
%PROJ4INV   Inverse map projection using PROJ.4 library
%   [LAT,LON]=PROJ4INV(PROJ,X,Y) returns the latitude and longitude
%   (degrees) from the inverse projection.  PROJ is a string with the PROJ.4 
%   +args syntax or a name of a supported projection. X and Y are arrays 
%   with map coordinates. 
%
%   To print supported projections use proj4defs();
%
%   Example:
%
%       [lat,lon]=proj4inv('+proj=utm +lon_0=112w +ellps=clrk66',460769.27,5011648.45)
%       result = 
%       45.25919447	-111.50000000 
%
%   See also proj4fwd, cs2cs and proj4defs.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015.

PROJEXE=proj4defs('PROJEXE');

if nargin~= 3
  error('Function must have three input arguments')
end
if length(x) ~= length(y)
  error('X and Y must have same length')
end

proj=proj4defs(proj);

infile=tempname;
fid=fopen(infile,'w');
for k=1:length(x)
   fprintf(fid,' %+16.4f  %+16.4f \n',x(k),y(k));
end
fclose(fid);
[status,result] = system( [ '"' PROJEXE '" -f "%.8f" -s -I ' proj ' ' infile ] );
delete(infile);
if status ~= 0
  error(result)
end
if nargout == 0
  disp(result)
elseif nargout == 1
  lat=result;
else
  c=textscan(result,'%f %f');
  lat=c{1};
  lon=c{2};
end

end


