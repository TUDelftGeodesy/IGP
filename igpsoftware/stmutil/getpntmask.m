function pntMask=getpntmask(pntCrd,roi) 
%GETPNTMASK  Compute point mask from Region Of Interest (ROI) file/string/numerical. 
%   PNTMASK=GETPNTMASK(PNTCRD,ROI) compute the point mask PNTMASK for the
%   coordinates in PNTCRD from the Region of Interest (ROI). PNTCRD is a 
%   column matrix with the coordinates (first two columns are used). The 
%   ROI can be 
%   - a string with the name of a shape or KML file containing the polygon.
%     The extension must be .shp or .kml.
%   - a 2-by-2 bounding box array
%   - a n-by-2 array with a polygon
%   The polygon, shapefile and kml file can contain multiple polygons with 
%   external and inner loops (holes). The external and internal loops should 
%   have opposite orientations; for example, a counterclockwise outer loop 
%   and clockwise inner loops or vice versa.
%
%   Examples:
%
%      pntMask = getpntmask(pntCrd,'example.shp')              % Shapefile with polygons
%      pntMask = getpntmask(pntCrd,'example.kml')              % KML file with polygons
%      pntMask = getpntmask(pntCrd,[53.25 6.5 ; 53.46 6.9 ])   % Bounding box  
%
%   See also ROI2POLY and INPOLYGON.
%
%   (c) Hans van der Marel, Delft University of Technology, 2021. 

% Created:  22 Oct 2021 by Hans van der Marel
% Modified: 
%

if isempty(roi)
   pntMask=true(size(pntCrd,1),1);
else
   poly=roi2poly(roi);
   pntMask=inpolygon(pntCrd(:,1),pntCrd(:,2),poly(:,1),poly(:,2)); 
end

end