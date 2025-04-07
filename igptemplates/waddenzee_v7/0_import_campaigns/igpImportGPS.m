%% Integrated geodetic processing - Import GPS levelling data
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project levelling data import module.
%
% Script to import levelling datasets and create space-time matrices.
%
% Input file
% - GPS levelling data in Cupido format 
%
% Ouputfile
% - space time dataset with GPS levelling data
%

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
igpimport('crsutil');
igpimport('rdnaptrans');

%% Set the input filenames, output filename and options

% Input file names and output file name

%inputfile={'../../../igpdata/nam_waddenzee/lts2_allgps_v6.0_20220111.nc'}; %format {'file1.nc';file2.nc';...}
%outputfile='lts2_allgps_v6.0_20220111.mat';
inputfile={'../../../igpdata/nam_waddenzee/lts2_allgps_v7.0_20230116.nc'};
outputfile='lts2_allgps_v7.0_20230116.mat';

% Options for Cupido conversion to stm

% Candidates for excluding campaigns
%
% prjName               dYear     yyyy-mmm-dd  dDays  #Pnts
% --------------------  --------  -----------  -----  -----
% NAM_GPS11P            2011.741  2011-Sep-29      6      5
% NAM_GPS15O            2015.517  2015-Jul-09      8     15
% NAM_GPS18O            2018.460  2018-Jun-18      5     15
%
% NAM_GPS08B            2008.801  2008-Oct-20     77      5
%
% NAM_GPS12L            2012.467  2012-Jun-20     30     11
% NAM_GPS15L            2015.495  2015-Jul-01     35     14
% NAM_GPS18L            2018.446  2018-Jun-13     28     15

options=[];
%options.POI=[2006,+Inf];
options.excludePrjName={ 'NAM_GPS11P' 'NAM_GPS15O' 'NAM_GPS18O' };
options.sdObsFlag = [0];                                    % Use only observations with sdobsflag 0
options.heightOnly=true;                                    % Only use the GPS up component!
%options.globalAttrib=globalAttribDefault;                  % Use defaults defined by igpinit
options.globalAttrib=[];                                    % Copy global attributes from Cupido NetCDF (is default)
options.globalAttrib.comment='Converted from Cupido/NetCDF';% Add comment attributes copied from NetCDF 

%% Do the conversion and save space time matrix to disk

cupido2stm(inputfile,outputfile,options,'overwrite');

%% What to do next ...?
%
%  1. Plot the network, timeseries and/or velocities
%  2. Integrate the 2006+ levelling data, GPS CORS and InSAR 
%
% Examples of plotting:
%
%   % cd to output folder
%   cd ../9_output
%   
%   % Plot network map by epoch -> plots are saved to a pdf 
%   stmplotmapbyepoch('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','saveplt','./');
%
%   % Estimate velocities and plot (using default 'refsystem'='min-velocities'
%   stmplotvel('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','saveplt','./');
%
%   % Plot the timeseries
%   stmplotseries('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','schranking','none','marker',true,'saveplt','./');
%   stmplotseries('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','schranking','ssdnorm','marker',true);
%
%   Notes: 
%
%   - plots are saved to pdf in the current folder, if you don't want this, leave out the 'saveplt' option, you can 
%     also specify another folder in 'saveplt'
%
%   - if you want to run from another directory than ../9_ouput, add options to find the some of the files: 
%
%     stmplotmapbyepoch(...,'background','../9_output/igp-background.mat')
%     stmplotvel(...,'background','../9_output/igp-background.mat','unstablearea','../9_ouput/igp-unstablearea.mat')
%     stmplotseries(...,'unstablearea','../9_ouput/igp-unstablearea.mat');
%
%   - you may wish to explore some other options, see the build in Matlab help
%
% [End of script]
