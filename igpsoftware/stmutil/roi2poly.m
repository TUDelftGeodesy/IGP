function output=roi2poly(roi,outformat)
%ROI2POLY  Convert Region Of Interest (ROI) file/string/numerical into a polygon. 
%   OUTPUT=ROI2POLY(ROI,OUTFORMAT) converts Region of Interest (ROI) file/string/numerical
%   into a Matlab polygon or geostruct. The ROI can be 
%   - a string with the name of a shape, KML or coordinate file containing the polygon.
%     The extension must be .shp, .kml or .coo.
%   - a 2-by-2 bounding box array
%   - a n-by-2 array with a polygon
%   The OUTPUT polygon can contain multiple polygons with external and inner 
%   loops (holes). The external and internal loops should have opposite 
%   orientations; for example, a counterclockwise outer loop and clockwise 
%   inner loops or vice versa. The output format OUTFORMAT can either be 'poly'
%   (default) or 'geostruct'.
%
%   Examples:
%
%      poly = roi2poly('example.shp')              % Shapefile with polygons
%      poly = roi2poly('example.kml','poly')       % KML file with polygons
%      poly = roi2poly([53.25 6.5 ; 53.46 6.9 ])   % Bounding box
%      poly = roi2poly('example2.coo','geostruct') % .coo file 
%
%   (c) Hans van der Marel, Delft University of Technology, 2021. 

% Created:  22 Oct 2021 by Hans van der Marel
% Modified: 12 Mar 2022 by Freek van Leijen
%              - enabled output of geostruct.
%           17 Jul 2023 by Freek van Leijen
%              - added reading of NAM .coo files
%            7 Mar 2024 by Hans van der Marel
%              - restored backwards compatibility


if nargin < 2 || isempty(outformat)
  outformat = 'poly';
end

if isnumeric(roi)

    if strmatch(outformat,'poly')
      output = bbox2poly(roi);
    elseif strmatch(outformat,'geostruct')
      output = bbox2geostruct(roi);  %%% to be implemented
    end
elseif ischar(roi) || istring(roi)

    % get geostruct
    geostruct = getroi(roi);

    if strmatch(outformat,'poly')
      if isfield(geostruct,'Lon')
        % check geostruct for polygons
        lat=[];
        lon=[];
        for k=1:numel(geostruct)
          if strcmpi(geostruct(k).Geometry,'Polygon')
             lat = [ lat ; geostruct(k).Lat(:) ];
             lon = [ lon ; geostruct(k).Lon(:) ];
          end
        end
        % convert to simple polygon
        output = [ lat lon ];
      elseif isfield(geostruct,'X')
        % check geostruct for polygons
        x=[];
        y=[];
        for k=1:numel(geostruct)
          if strcmpi(geostruct(k).Geometry,'Polygon')
             x = [ x ; geostruct(k).X(:) ];
             y = [ y ; geostruct(k).Y(:) ];
          end
        end
        % convert to simple polygon
        output = [ x y ];
      else
        error('The geostruct has a wrong format.');
      end

    elseif strmatch(outformat,'geostruct')
      output = geostruct;

    else
      error('You specified a wrong ROI format.');
    end

else
   error('roi must be numeric or character type')
end
 

end


function poly=bbox2poly(bbox)
%BBOX2POLY  Convert bounding box into a polygon.
%   POLY=BBOX2POLY(BBOX) converts bounding box BBOX into a polygon POLY.
%   BBOX can be a bounding box or polygon. The function checks for proper
%   closure of the polygon. POLY can be used by Matlab's INPOLYGON().
%
%   Examples:
%
%      poly = bbox2poly([53.25 6.5 ; 53.46 6.9 ])
%
%   (c) Hans van der Marel, Delft University of Technology, 2020. 

% Created:  8 Oct 2020 by Hans van der Marel
% Modified: 

if size(bbox,1) == 2 && size(bbox,2) == 2 
   % make a polygon out of bounding box
   poly=[ bbox(1,1) bbox(1,2) ; ...    
          bbox(2,1) bbox(1,2) ; ...
          bbox(2,1) bbox(2,2) ; ...
          bbox(1,1) bbox(2,2) ; ...
          bbox(1,1) bbox(1,2) ];    
elseif size(bbox,1) > 2 && size(bbox,2) == 2
   % bbox is already a polygon
   poly=bbox;    
else
   error('This function expects a 2x2 bounding box  [latmin lonmin; latmax lonmax] or polygon.')    
end

% Close the polygon
if poly(1,1) ~= poly(end,1) || poly(1,2) ~= poly(end,2)
   poly = [poly ; poly(1,:) ];
end

end

function roi = getroi(roiFile)
%GETROI Get the Region Of Interest (ROI) from file. 
%   GETROI(ROIFILE) Get the Region of Interest (ROI) based on
%   a .shp or .kml file containing a polygon. The output is
%   a geostruct.
%
%   roi = getroi(roiFile)
%
%   Examples:
%
%      roi = getroi('example.shp');
%      roi = getroi('example.kml');
%
%   TODO: implement .json, .wkt formats
%
%   (c) Freek van Leijen, Delft University of Technology, 2021. 

% Created:  12 Jan 2021 by Freek van Leijen
% Modified: 22 Oct 2021 by Hans van der Marel
%              - use as internal function for ROI2POLY

[~,~,extension] = fileparts(roiFile);

switch extension
  case '.shp'
    roi = shaperead(roiFile, 'UseGeoCoords', true);
  case '.kml'
    roi = kml2struct(roiFile);
    %roi = kml2struct_multi(roiFile); % Would be needed for MULTIPOLGONS, but does not work properly
  case '.coo'
    roi = read_polygon(roiFile);
  case {'.json','.geojson'}
    error('The JSON format is not implemented yet.');
    %roi = jsondecode(fileread(roiFile));  %TODO reshape struct
  case {'.wkt'}
    error('The WKT format is not implemented yet.');
    %% Create WKT at https://arthur-e.github.io/Wicket/sandbox-gmaps3.html
    %fid = fopen(roiFile,'r');
    %wkts = textscan(fid,'%s','delimiter','#'); % # is dummy value
    %fclose(fid);
    %roi = wkt2geostruct(wkts{1}); %TODO check why order of points is reversed
  otherwise
    error('Your ROI file has an unsupported format/extension.');
end

end

function kmlStruct = kml2struct(kmlFile)
% kmlStruct = kml2struct(kmlFile)
%
% Import a .kml file as a vector array of shapefile structs, with Geometry, Name,
% Description, Lon, Lat, and BoundaryBox fields.  Structs may contain a mix
% of points, lines, and polygons.
%
% .kml files with folder structure will not be presented as such, but will
% appear as a single vector array of structs.
%
% 

[FID msg] = fopen(kmlFile,'rt');

if FID<0
    error(msg)
end

txt = fread(FID,'uint8=>char')';
fclose(FID);

expr = '<Placemark.+?>.+?</Placemark>';

objectStrings = regexp(txt,expr,'match');

Nos = length(objectStrings);

for ii = 1:Nos
    % Find Object Name Field
    bucket = regexp(objectStrings{ii},'<name.*?>.+?</name>','match');
    if isempty(bucket)
        name = 'undefined';
    else
        % Clip off flags
        name = regexprep(bucket{1},'<name.*?>\s*','');
        name = regexprep(name,'\s*</name>','');
    end
    
    % Find Object Description Field
    bucket = regexp(objectStrings{ii},'<description.*?>.+?</description>','match');
    if isempty(bucket)
        desc = '';
    else
        % Clip off flags
        desc = regexprep(bucket{1},'<description.*?>\s*','');
        desc = regexprep(desc,'\s*</description>','');
    end
    
    geom = 0;
    % Identify Object Type
    if ~isempty(regexp(objectStrings{ii},'<Point', 'once'))
        geom = 1;
    elseif ~isempty(regexp(objectStrings{ii},'<LineString', 'once'))
        geom = 2;
    elseif ~isempty(regexp(objectStrings{ii},'<Polygon', 'once'))
        geom = 3;
    end
    
    switch geom
        case 1
            geometry = 'Point';
        case 2
            geometry = 'Line';
        case 3
            geometry = 'Polygon';
        otherwise
            geometry = '';
    end
    
    % Find Coordinate Field
    bucket = regexp(objectStrings{ii},'<coordinates.*?>.+?</coordinates>','match');
    % Clip off flags
    coordStr = regexprep(bucket{1},'<coordinates.*?>(\s+)*','');
    coordStr = regexprep(coordStr,'(\s+)*</coordinates>','');
    % Split coordinate string by commas or white spaces, and convert string
    % to doubles
    coordMat = str2double(regexp(coordStr,'[,\s]+','split'));
    % Rearrange coordinates to form an x-by-3 matrix
    [m,n] = size(coordMat);
    coordMat = reshape(coordMat,3,m*n/3)';
    
    % define polygon in clockwise direction, and terminate
    [Lat, Lon] = poly2ccw(coordMat(:,2),coordMat(:,1));
    if geom==3
        Lon = [Lon;NaN];
        Lat = [Lat;NaN];
    end
    
    % Create structure
    kmlStruct(ii).Geometry = geometry;
    kmlStruct(ii).Name = name;
    kmlStruct(ii).Description = desc;
    kmlStruct(ii).Lon = Lon;
    kmlStruct(ii).Lat = Lat;
    kmlStruct(ii).BoundingBox = [[min(Lon) min(Lat);max(Lon) max(Lat)]];
end

end

function out = read_polygon(filename)
%% READ_POLYGON(filename)
%   usage: xy = read_polygon(filename,...)
%
%   where filename indicates a file in the following format:
%
%   x11 y11  % first disconnected polygon
%   x12 y12
%   x13 y13
%   x11 y11  % repeat first node
%
%   x21 y21  % second disconnected polygon
%   x22 y22
%   ...

fid = fopen(filename);
if fid<0
    fprintf(1,'File %s not found.\n',filename);
end

data = textscan(fid,'%s','delimiter','\n');
fclose(fid);

breaks = [0;find(ismember(data{1},''))];
n = length(breaks);
breaks = [breaks;length(data{1})+1];
outlines = cell(n,1);

geometry = 'Polygon';
for ii = 1:n

  outlines{ii,1} = cellfun(@(s) textscan(s,'%f %f','delimiter',' ,'),...
       data{1}(breaks(ii)+1:breaks(ii+1)-1),'uniformoutput',false);
  outlines{ii,1} = cell2mat(vertcat(outlines{ii}{:}));
   
  % Create structure
  X = outlines{ii}(:,1)';
  Y = outlines{ii}(:,2)';
  out(ii).Geometry = geometry;
  out(ii).X = X;
  out(ii).Y = Y;
  out(ii).BoundingBox = [min(X) min(Y);max(X) max(Y)];

end


end

