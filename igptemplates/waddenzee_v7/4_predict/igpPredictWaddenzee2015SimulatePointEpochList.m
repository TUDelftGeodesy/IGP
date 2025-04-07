%% Integrated geodetic processing - Simulation of random Point and Epoch list
%% (for testing)
%
% *Freek van Leijen*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)
% project simulation of random Point and Epoch list (for testing).

%% Add required toolboxes to Matlab path
%
% 1. run script in 'igpproject' ('igpdata') project root to set Matlab path 
%    to  'igpsoftware' (so that we can find the function 'igpimport')
% 2. add the required toolboxes to the Matlab path using 'igpimport'
%    function
%
% 'igpinit' is a script that resides in the igpproject/igpdata  root, 'igpimport'
% is a function that resides in the 'igpsoftware' software directory. The
% location of the toolboxes is defined in 'igptoolbox.cfg' that is located
% in the same directory as 'igpimport'. It is possible to use alternative
% environments (e.g. for development) by specifying a new configuration
% file with paths to the toolboxes as second argument of 'igpimport'.

run ../../igpinit         % add igpsoftware folder to Matlab path 

igpimport('stmmain');     % add all required toolboxes to the Matlab path
igpimport('stmutil');
igpimport('proj');
igpimport('crsutil');
igpimport('rdnaptrans')

%inputFilename = '../simtest2a_combined_truth.mat';
inputFilename = '../3_integrate/waddenzee_2015.mat';

%crstoolbox = 'proj4';
crstoolbox = 'rdnaptrans';
crsIn= 'wgs84';
crsOut = 'rd';

Npoints = 40;
Nepochs = 24;

stmIn = stmread(inputFilename);

latMin = min(stmIn.pntCrd(:,1));
latMax = max(stmIn.pntCrd(:,1));
lonMin = min(stmIn.pntCrd(:,2));
lonMax = max(stmIn.pntCrd(:,2));

lat = rand(Npoints,1)*(latMax-latMin)+latMin;
lon = rand(Npoints,1)*(lonMax-lonMin)+lonMin;

fid = fopen('pointlist_wgs84_waddenzee_2015.txt','w');
fprintf(fid,['%12.8f %12.8f\n'],[lat lon]');
fclose(fid);

lat2 = rand(Npoints,1)*(latMax-latMin)+latMin;
lon2 = rand(Npoints,1)*(lonMax-lonMin)+lonMin;

rd = crstrans([lat lon],crsIn,crsOut,crstoolbox);

fid2 = fopen('pointlist_rd_waddenzee_2015.txt','w');
fprintf(fid2,['%12.3f %12.3f\n'],rd(:,1:2)');
fclose(fid2);

epochMin = min(stmIn.epochDyear);
epochMax = max(stmIn.epochDyear);

epochs1 = rand(Nepochs,1)*(epochMax-epochMin)+epochMin;
epochs2 = rand(Nepochs,1)*(epochMax-epochMin)+epochMin;

fid3 = fopen('epochlist1_waddenzee_2015.txt','w');
fprintf(fid3,['%6.2f\n'],epochs1);
fclose(fid3);

fid4 = fopen('epochlist2_waddenzee_2015.txt','w');
fprintf(fid4,['%6.2f\n'],epochs2);
fclose(fid4);

