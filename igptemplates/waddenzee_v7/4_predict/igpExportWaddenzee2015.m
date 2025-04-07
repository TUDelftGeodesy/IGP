%% Integrated geodetic processing - Waddenzee export
%
% *Freek van Leijen*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)
% project Waddenzee export.

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

dbstop if error

%% Export 
%
% Export of results for NAM Integrated Geodetic Processing (IGP) project
%
% Syntax
% 
%    stmexport(inputFilename,options)
%
% Options can be specified as a structure 
%
%   options=[];
%   options.format={ 'geotiff' 'csv' 'png' 'shp' };
%   options.unit = 'm';
%
% or cell array
%
%    options={'format',{'csv' 'png'},,'unit','m'};
%


%% Set STM file
options = [];

%% Set output format (can be multiple)
%options.format = { 'csv' 'shp' };
%options.format = { 'png' 'geotiff' };
options.format = { 'png' 'csv' };
%options.format = { 'geotiff' }:
%options.format = { 'geotiff' 'csv' 'png' 'shp' };
%options.format = { 'csv' };

%% Set plots
options.dpi = '300';
options.cmap = 'defo';
options.visible = 'on';
options.background = '../9_output/igp-background-rd.mat';

%% Set output unit
%options.unit = 'mm';
options.unit = 'm';

%% Perform export
inputFilename = './waddenzee_2015_predicted.mat';
stmexport(inputFilename,options);

options.background = '../9_output/igp-background.mat';
inputFilename = './waddenzee_2015_predicted_list_wgs84.mat';
stmexport(inputFilename,options);

options.background = '../9_output/igp-background-rd.mat';
inputFilename = './waddenzee_2015_predicted_list_rd.mat';
stmexport(inputFilename,options);

inputFilename = './waddenzee_2015_predicted_list_rd_roi.mat';
stmexport(inputFilename,options);

% [End of script]

