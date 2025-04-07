%% IGP contour of unstable area
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)

%% Add required toolboxes to Matlab path

run ../../igpinit         % add igpsoftware folder to Matlab path 

igpimport('crsutil');     % add all required toolboxes to the Matlab path
igpimport('rdnaptrans');

%% Read polygon defining unstable area from .coo text file with RD coordinates

coofile='unstablearea.coo';

if exist('readmatrix','file')
   unstable_area_rd = readmatrix(coofile,'FileType','text','emptyLineRule','read');
else
   unstable_area_rd = [];
   fid=fopen(coofile);
   while ~feof(fid)
      tline = fgetl(fid);
      crd = sscanf(tline,'%f');
      if numel(crd) ~= 2
          crd = nan(1,2);
      end
      unstable_area_rd=[unstable_area_rd ; crd(:)'];
   end
   fclose(fid);
end

%% Convert to latitude / longitude matrix

unstable_area = rd2etrs(unstable_area_rd);  % <- required modified rd2etrs, needs option 'ZERO'

%% Plot the resulting data (with igp-background)

figure
if exist('igp-background.mat','file')
   igpbackground=load('igp-background');
   plot(igpbackground.Lon,igpbackground.Lat)
   hold on
end
plot(unstable_area(:,2),unstable_area(:,1))
grid on
daspect([1/cosd(nanmean(unstable_area(:,1))) 1 1])
xlabel('Longitude [deg]')
ylabel('Latitude [deg]')


%% Create a matfile igp-unstablearea.mat with lat lon data

Lat=unstable_area(:,1);
Lon=unstable_area(:,2);

save('igp-unstablearea','Lat','Lon')

disp('Move igp-unstablearea.mat to 9_output or a directory in the search path')

