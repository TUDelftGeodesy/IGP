function [x,y]=proj4fwd(proj,lat,lon)
%PROJ4FWD   Forward map projection using PROJ.4 library
%   [X,Y]=PROJ4FWD(PROJ,LAT,LON) returns the X and Y map coordinates
%   from the forward projection.  PROJ is a string with the PROJ.4 
%   +args syntax or a name of a supported projection. LAT and LON are 
%   arrays with latitude and longitude in degrees. 
%
%   To print supported projections use proj4defs();
%
%   Example:
%
%       [x,y]=proj4fwd('+proj=utm +lon_0=112w +ellps=clrk66',+45.25919447,-111.5)
%       result = 
%       460769.27   5011648.45
%
%   See also proj4inv, cs2cs and proj4defs.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015.

PROJEXE=proj4defs('PROJEXE');

if nargin~= 3
  error('Function must have three input arguments')
end
if length(lat) ~= length(lon)
  error('LAT and LON must have same length')
end

proj=proj4defs(proj);

infile=tempname;
fid=fopen(infile,'w');
for k=1:length(lat)
   fprintf(fid,' %+16.11f  %+16.11f \n',lon(k),lat(k));
end
fclose(fid);
[status,result] = system( [ '"' PROJEXE '" ' proj ' ' infile ] );
delete(infile);
if status ~= 0
  error(result)
end
if nargout == 0
  disp(result)
elseif nargout == 1
  x=result;
else
  c=textscan(result,'%f %f');
  x=c{1};
  y=c{2};
end

end


