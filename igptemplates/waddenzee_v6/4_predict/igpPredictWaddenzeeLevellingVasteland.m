%% Integrated geodetic processing - Waddenzee prediction
%
% *Freek van Leijen*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)
% project Waddenzee prediction.

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

dbstop if error

%% Prediction 
%
% Prediction for NAM Integrated Geodetic Processing (IGP) project
%
% Syntax
% 
%    igppredict(inputFilename,outputFilename,options,flag)
%
% Options can be specified as a structure 
%
%   options=[];
%   options.verbose=1;
%   options.crstoolbox = 'PROJ4';
%
% or cell array
%
%    options={'verbose',1,'crstoolbox','PROJ4'};
%
%   flag='overwrite';
%


%% Set STM file
inputFilename = '../0_import_campaigns/lts2_levelling_vasteland_full_period.mat';
%outputFilename = './waddenzee_2015_predicted.mat'; %Various output files defined below
options = [];

%% Set coordinate transformation toolbox
options.crstoolbox = 'RDNAPTRANS';
%options.crstoolbox = 'rdnaptrans';
%options.crstoolbox = 'proj4';
%options.crstoolbox = 'PROJ4';

%% Set Region of Interest
%options.ROI = 'area.shp'; % .kml, .shp
%options.ROI = 'groningen_test_roi.kml';
%options.ROI = 'pointlist.txt'; % .txt, .csv, .xlsx
%options.ROI = [52.3434 6.34343; 52.2534 6.34342];
%options.ROI = [230000 570000; 250000 590000];
options.ROI = [];

%% Set CRS ROI
%options.crsROI = 'UTM';
%options.crsROI = 'WGS84';
%options.crsROI = 'RD';
%options.crsROI = 'EPSG28992';
options.crsROI = [];

%% Set Period of Interest (in dates or decimal years)
% periods {{start,stop},{start,stop}, ...}, epoch list (cellarray,.txt, .csv, .xlsx)
%options.POI = { {'1984-01-01','2021-01-01'} }; % Date format should be 'YYYY-MM-DD' (for now)
%options.POI = { {'1984-01-01','1990-01-01'} {'2000-01-01','2021-01-01'} };
%options.POI = [ 1984.00 2021.00 ];
%options.POI = [ 1984.00 1990.00; 2000.00 2021.00 ];
options.POI = [];

%% Set Spatial Reference
%options.spatialRef = [];
options.spatialRef = {'../9_output/unstablearea.coo' 'EPSG28992' 'reverse'}; %
%options.spatialRef = {'filename.kml' 'wgs84' 'normal'};

%% Set Temporal Reference
options.temporalRef = [];
%options.temporalRef = 1984;

%% Set CRS output
%options.crsOutput = 'UTM';
%options.crsOutput = 'WGS84';
options.crsOutput = 'RD';
%options.crsOutput = 'EPSG28992';

%% Set spatial output
%options.spaceOutput = { 'original' };
options.spaceOutput = { 'grid' 500 'm' };
%options.spaceOutput = { 'grid' 0.001 'deg' };
%options.spaceOutput = { 'list' 'pointlist.txt' 'm' };

%% Set temporal output
%options.timeOutput = { 'original' };
options.timeOutput = { 'sampling' 1 'year' }; % e.g., 0.2 year
%options.timeOutput = { 'list' 'epochlist.txt' 'year' };

%% Set trend estimation settings
options.ignoreStochModel = false;
options.refsystem = 'min-velocities';

%% Set covariance function of signal ('estimated' or one model per obs type)
% Use double cells, can be multiple models
% e.g., options.signalModel = { {'[modelstring]' '[modelstring]'} {'[modelstring]'} };
%options.signalModel = { 'estimated' }; %NOT IMPLEMENTED YET
%options.signalModel = {{'tudpred1(s20=NaN,s2t=3,s2s=3,Rt=1,Rs=5)'} ...
%                       {'tudpred1(s20=NaN,s2t=2,s2s=2,Rt=1,Rs=3)'}};
%options.signalModel = {{'tudpred1(s20=NaN,s2t=4,s2s=3,Rt=1,Rs=4)'} ...
%                       {'tudpred1(s20=NaN,s2t=2,s2s=2,Rt=1,Rs=3)'}};
%options.signalModel = {{'sure(s2z=0.8,L=4,q=0.85,tref=2003)'} ...
%                       {'sure(s2z=0.8,L=4,q=0.85,tref=2003)'} ...
%                       {'sure(s2z=0.8,L=4,q=0.85,tref=2003)'} ...
%                       };
%options.signalModel = {{'sure(s2z=0.8,L=4,q=0.85,tlag=0)'} ... % tlag=1 means tref becomes 1 year earlier, hence 1984 -> 1983
%                       {'sure(s2z=0.8,L=4,q=0.85,tlag=0)'} ...
%                       {'sure(s2z=0.8,L=4,q=0.85,tlag=0)'} ...
%                       };
options.signalModel = {{'sure(s2z=0.8,L=4,q=0.85,tlag=0)'}}; ... % tlag=1 means tref becomes 1 year earlier, hence 1984 -> 1983
%options.signalModel = {{'sure(s2z=0.8,L=4.0,q=0.85,tref=2010)'} ...
%                       {'houtenbos(s2t=0.16,pt=1.24)'}}; 

%% Set clustering settings
options.clusteringMethod = 'distance';
options.spatialTolerance = 0.1;
options.doplots = 0;


%% Set prediction settings
options.predConvThres = 1000;
options.predMaxIter = 1;

%% Set interpolation settings (for trend interpolation)
%options.interpMethod = 'linear';
options.interpMethod = 'idw';
options.idwPower = 2;
%options.idwPower = 4;

%% Set plot settings
options.background = '../9_output/igp-background.mat';
options.plotmode = 'last';
%options.saveplt = [];
options.saveplt = 'generated';

%% Perform prediction

outputFilename = './lts2_levelling_vasteland_full_period_predicted.mat';
stmpredict(inputFilename,outputFilename,options,'overwrite');

outputFilename = './lts2_levelling_vasteland_full_period_predicted_list_wgs84.mat';
options.crsOutput = 'wgs84';
options.spaceOutput = { 'list' 'pointlist_wgs84_levelling_vasteland.txt' 'deg' };
options.timeOutput = { 'list' 'epochlist1_levelling_vasteland.txt' 'year' };
stmpredict(inputFilename,outputFilename,options,'overwrite');

outputFilename = './lts2_levelling_vasteland_full_period_predicted_list_rd.mat';
options.crsOutput = 'RD';
options.spaceOutput = { 'list' 'pointlist_rd_levelling_vasteland.txt' 'm' };
options.timeOutput = { 'list' 'epochlist2_levelling_vasteland.txt' 'year'};
stmpredict(inputFilename,outputFilename,options,'overwrite');

outputFilename = './lts2_levelling_vasteland_full_period_predicted_list_rd_roi.mat';
options.ROI = [230000 570000; 250000 590000];
options.crsROI = 'EPSG28992';
%options.crstoolbox = 'proj4';
stmpredict(inputFilename,outputFilename,options,'overwrite');

% [End of script]

