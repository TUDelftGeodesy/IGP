%% IGP background map 
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)

%% Add required toolboxes to Matlab path

run ../../igpinit         % add igpsoftware folder to Matlab path 

igpimport('basemaps');    % add all required toolboxes to the Matlab path
igpimport('rdnaptrans');

%% Select region of interest and plot using mapping toolbox

latrange=[53.0 53.6];
lonrange=[5.5 7.25 ];

nweuropemap(latrange,lonrange,'F')

%% Prepare matfile with Lon/Lat polygons 

% read shapefiles (as geostruct files)

id='nweurope';
resolution='f';

bbox = [ lonrange' latrange'];

landareas=shaperead(['landareas_' id '_' resolution],'UseGeoCoords', true, 'BoundingBox', bbox);
lakes=shaperead(['lakes_' id '_' resolution],'UseGeoCoords', true, 'BoundingBox', bbox);
islands_in_lakes=shaperead(['islands_' id '_' resolution],'UseGeoCoords', true, 'BoundingBox', bbox);
borders=shaperead(['borders_' id '_' resolution],'UseGeoCoords', true, 'BoundingBox', bbox);
rivers=shaperead(['rivers_' id '_' resolution],'UseGeoCoords', true, 'BoundingBox', bbox);
 
% Convert geostruct to line segments and crop the data

multiline=geostruct2line(landareas,lakes,borders);
multiline=cropmultiline(multiline,latrange,lonrange);

% Plot the resulting data

figure
plot(multiline(:,2),multiline(:,1))
xlim(lonrange)
ylim(latrange)
grid on
daspect([1/cosd(mean(latrange)) 1 1])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')


%% Create a matfile igp-background.mat with Lat and Lon arrays

Lat=multiline(:,1);
Lon=multiline(:,2);

%save('igp-background','Lat','Lon')

disp('Move igp-background.mat to 9_output or a directory in the search path')

nanIdx = ~isnan(Lat);
XY = crstrans([Lat(nanIdx) Lon(nanIdx)],'wgs84','rd','rdnaptrans');
X = NaN(size(Lat));
Y = NaN(size(Lat));
X(nanIdx) = XY(:,1);
Y(nanIdx) = XY(:,2);

save('igp-background-rd','X','Y')

disp('Move igp-background-rd.mat to 9_output or a directory in the search path')

%% Local functions

function multiline=geostruct2line(varargin)

multiline=[];
for l=1:nargin
   geostruct=varargin{l};
   for k=1:numel(geostruct)
      multiline = [ multiline ; geostruct(k).Lat' geostruct(k).Lon' ; nan nan];
   end
end

end

function multiline=cropmultiline(multiline,latrange,lonrange)

mask = multiline(:,1) >= latrange(1) & multiline(:,1) <= latrange(2) &  ...
       multiline(:,2) >= lonrange(1) & multiline(:,2) <= lonrange(2) ;
multiline(~mask,:)=nan;

nanrows=all(isnan(multiline),2);
nanrows= nanrows & [ diff(nanrows) ; 0 ] == 0;

multiline(nanrows,:)=[];

end

