function [xout,yout]=cs2cs(projin,projout,xin,yin)
%CS2CS   Coordinate transformation using PROJ.4 library
%   [XOUT,YOUT]=CS2CS(PROJIN,PROJOUT,XIN,YIN) returns map coordinates 
%   XOUT and YOUT in the cartographic coordinate projection PROJOUT.
%   XIN and YIN are arrays with map coordinates in the projection PROJIN.
%   PROJIN and PROJOUT are strings with the PROJ.4 +args syntax and/or 
%   names of supported projections. 
%
%   [XOUT,YOUT]=CS2CS(PROJIN,PROJOUT,INFILE) reads the input coordinates
%   from the file INFILE.
%
%   RESULT=CS2CS(...) return the result as text.
%
%   CS2CS(...) print the result on standard outreturn the result as text.
%
%   CS2CS performs transformation between the source and destination 
%   cartographic coordinate system on a set of input points. The coordinate
%   system transformation can include translation between projected and 
%   geographic coordinates as well as the application of datum shifts.
%
%   In case input and/or output coordinates are latitude and longitude
%   X corresponds to latitude in degrees and Y to longitude in degrees.
%
%   To print supported projections use proj4defs();
%
%   Example:
%
%       [x,y]=cs2cs('+proj=latlong +datum=NAD83','+proj=utm +zone=10 +datum=NAD27',+45.25919447,-111.5)
%       result =
%       1402224.58	5076275.43 -0.00
%       [x,y]=cs2cs('+proj=utm +zone=10 +datum=NAD27','+proj=latlong +datum=NAD83',1402224.58,5076275.43)
%       result =
%       45.25919435	-111.50000123 0.00005280 
%
%       will  transform the input NAD83 geographic coordinates into NAD27 
%       coordinates in the UTM projection with zone 10 selected, and vice
%       versa. 
% 
%   See also proj4fwd, proj4inv and proj4defs.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015.
%
%   Modified: 18 Jan 2021 by Freek van Leijen
%             - added formatting for '+proj=longlat' output projection

CS2CSEXE=proj4defs('CS2CSEXE');

if nargin < 3 || nargin > 4
  error('Function must have three or four input arguments')
end
if nargin == 4
  if length(xin) ~= length(yin)
    error('XIN and YIN must have same length')
  end
end

projin=proj4defs(projin);
projout=proj4defs(projout);

options='';
if strfind(projin,'+proj=latlong')
  options=[ options '-r '];
end
if strfind(projout,'+proj=latlong')
  options=[ options '-f "%.8f" -s '];
elseif strfind(projout,'+proj=longlat')
  options=[ options '-f "%.8f" '];
end

if nargin == 4
  infile=tempname;
  fid=fopen(infile,'w');
  for k=1:length(xin)
     fprintf(fid,' %+20.11f  %+20.11f \n',xin(k),yin(k));
  end
  fclose(fid);
else
  infile=xin;
end
[status,result] = system( [ '"' CS2CSEXE '" ' options ' ' projin ' +to ' projout ' ' infile ] );
if nargin == 4
  delete(infile);
end
if status ~= 0
  error(result)
end

if nargout == 0
  disp(result)
elseif nargout == 1
  xout=result;
else
  c=textscan(result,'%f %f %f');
  xout=c{1};
  yout=c{2};
end

end


