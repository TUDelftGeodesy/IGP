%% Integrated geodetic processing - Import GPS levelling data
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project levelling data import module.
%
% Script to import levelling datasets and create space-time matrices.
%
% Input file
% - levelling data in Cupido format
%
% Ouputfiles 
% - space time datasets with subsets of levelling data:
%      - lts2_levelling_vasteland.mat
%      - lts2_levelling_vasteland_2006.mat
%      - lts2_levelling_ameland.mat
%      - lts2_levelling_ameland_2006.mat
%   Other subsets can be selected by uncommenting parts of the code
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

% Input file waddenzee

%inputfile={'../../../igpdata/nam_waddenzee/lts2_alllevelling_gpsbaseline_v5.2_20210323.nc'}; 
inputfile={'../../../igpdata/nam_waddenzee/lts2_alllevelling_gpsbaseline_v7.1_20230222.nc'};

% OVERVIEW OF OPTIONS FOR SELECTING NETWORKS
% ------------------------------------------
%
% See "Overview_of_levellingg_campaigns.xlsx" for a complete overview

% All campaigns

% options.includePrjName = { };
% options.excludePrjName={ };
% outputfile='lts2_levelling_all.mat';

% Hydrostatisch ?
%
% 279H05                1986.476  1986-Jun-24   1819     48  AM+2VL
% 344W01                1993.975  1993-Dec-23    178     31  AM+SC+2VL
% 371W00                1998.329  1998-May-01     91    195  AM+SC+xxVL
% 373W31                2001.504  2001-Jul-04    299     59  AM+SC+2VL
% options.includePrjName = { ... 
%   '279H05'  '344W01'  '371W00'  '373W31' ...
% };
% options.excludePrjName={ };
% outputfile='lts2_levelling_other.mat';

% GPS 
% options.includePrjName = { ... 
%  'NAM_GPS93' 'NAM_GPS97' 'NAM_GPS00' 'NAM_GPS04A' ...
% };
% options.excludePrjName={ };
% outputfile='lts2_levelling_bsl.mat';

% Ameland  (See "Overview_of_levellingg_campaigns.xlsx")
%
% options.includePrjName = { ... 
%   '289W05'  '289W16'  '289W20'  '289W26'  '289W34'  '279W22'  '289W37'  '332W02'  '342W04'  '342W05'  '342W10' ... 
%   'aml1998'  'aml1999' 'aml2003'  'AMEL0409' 'aml2005' ...
%  '70604-001' '231-70604' '231-70812' 'Amel_2009' 'NAM_AM2011' 'NAM_AM2014' 'NAM_AM2017' 'NAM_AM2020' ...
% };
% options.excludePrjName={ };
% outputfile='lts2_levelling_ameland.mat';

% Schiermonnikoog  (See "Overview_of_levellingg_campaigns.xlsx")
%
% options.includePrjName = { ... 
%   'NAM_SC2006'  'NAM_SC2009'  'NAM_SC2012'  'NAM_SC2015'  'NAM_SC2018'  'NAM_SC2021'  ...
% };
% options.excludePrjName={ };
% outputfile='lts2_levelling_schier.mat';

% Vasteland (See "Overview_of_levellingg_campaigns.xlsx")
%
% options.includePrjName = { ...
%   '231'  '300'  '289W22'  '289W31'  '332W04'  '342W08'  '342W12'  '364W00'  '365W37-SEC'  '370W27'  'GRON0306' ...
%   'NAM_GF2008' 'NAM_LM2011' 'NAM_LM2015'  'NAM_MA2006' 'RWS_GF2013'  'RWS_GF2018' ...
% };
% options.excludePrjName={ '231-'};
% outputfile='lts2_levelling_vasteland.mat';

% Options for ROI
% options.POI=[2006,+Inf];


%% All - full period

outputfile='lts2_levelling_v7_all.mat';

options=[];
options.splitNetwork='delete';
options.includePrjName = { };
options.excludePrjName={ };

options.sdObsFlag = [0];                                    % Use only observations with sdobsflag 0
%options.globalAttrib=globalAttribDefault;                  % Use defaults defined by igpinit
options.globalAttrib=[];                                    % Copy global attributes from Cupido NetCDF (is default)
options.globalAttrib.comment='Converted from Cupido/NetCDF';% Add comment attributes copied from NetCDF 

% Do the conversion and save space time matrix to disk

cupido2stm(inputfile,outputfile,options,'overwrite');


%% Vasteland 

outputfile='lts2_levelling_vasteland_v7_full_period.mat';

options=[];
options.splitNetwork='delete';
options.includePrjName = { ...
  '231'  '300'  '289W22'  '289W31'  '332W04'  '342W08'  '342W12'  '364W00'  '365W37-SEC'  '370W27'  'GRON0306' ...
  'NAM_GF2008' 'NAM_LM2011' 'NAM_LM2015'  'NAM_MA2006' 'RWS_GF2013'  'RWS_GF2018' ...
};
options.excludePrjName={ '231-'};

options.sdObsFlag = [0];                                    % Use only observations with sdobsflag 0
%options.globalAttrib=globalAttribDefault;                  % Use defaults defined by igpinit
options.globalAttrib=[];                                    % Copy global attributes from Cupido NetCDF (is default)
options.globalAttrib.comment='Converted from Cupido/NetCDF';% Add comment attributes copied from NetCDF 

% Do the conversion and save space time matrix to disk

cupido2stm(inputfile,outputfile,options,'overwrite');

%% Vasteland - from 2006 onwards

outputfile='lts2_levelling_vasteland.mat';
options.POI=[2006,+Inf];

cupido2stm(inputfile,outputfile,options,'overwrite');

%% Ameland 

outputfile='lts2_levelling_ameland_v7_full_period.mat';

options=[];
options.splitNetwork='delete';
options.includePrjName = { ... 
  '289W05'  '289W16'  '289W20'  '289W26'  '289W34'  '279W22'  '289W37'  '332W02'  '342W04'  '342W05'  '342W10' ... 
  'aml1998'  'aml1999' 'aml2003'  'AMEL0409' 'aml2005' ...
 '70604-001' '231-70604' '231-70812' 'Amel_2009' 'NAM_AM2011' 'NAM_AM2014' 'NAM_AM2017' 'NAM_AM2020' ...
};
options.excludePrjName={ };

options.sdObsFlag = [0];                                    % Use only observations with sdobsflag 0
%options.globalAttrib=globalAttribDefault;                  % Use defaults defined by igpinit
options.globalAttrib=[];                                    % Copy global attributes from Cupido NetCDF (is default)
options.globalAttrib.comment='Converted from Cupido/NetCDF';% Add comment attributes copied from NetCDF 

% Do the conversion and save space time matrix to disk

cupido2stm(inputfile,outputfile,options,'overwrite');

%% Ameland - from 2006 onwards

outputfile='lts2_levelling_ameland.mat';
options.POI=[2006,+Inf];

cupido2stm(inputfile,outputfile,options,'overwrite');

%% Schiermonnikoog

outputfile='lts2_levelling_schier_v7_full_period.mat';

options=[];
options.splitNetwork='delete';
options.includePrjName = { ... 
   'NAM_SC2006'  'NAM_SC2009'  'NAM_SC2012'  'NAM_SC2015'  'NAM_SC2018'  'NAM_SC2021'  ...
};
options.excludePrjName={ };

options.sdObsFlag = [0];                                    % Use only observations with sdobsflag 0
%options.globalAttrib=globalAttribDefault;                  % Use defaults defined by igpinit
options.globalAttrib=[];                                    % Copy global attributes from Cupido NetCDF (is default)
options.globalAttrib.comment='Converted from Cupido/NetCDF';% Add comment attributes copied from NetCDF 

% Do the conversion and save space time matrix to disk

cupido2stm(inputfile,outputfile,options,'overwrite');

%% Schiermonnikoog - from 2006 onwards

outputfile='lts2_levelling_schier.mat';
options.POI=[2006,+Inf];

cupido2stm(inputfile,outputfile,options,'overwrite');

%% What to do next ...?
%
%  1. Inspect stm, plot the levelling networks, timeseries and/or velocities
%  2. Use the prediction module on the full length datasets (once for each 
%     subnet)
%  3. Integrate the 2006+ data with GPS CORS and InSAR 
%
% Inspect the stm:
%
%   stmdisp('lts2_levelling_ameland.mat');
%   stmdisp('lts2_levelling_schier.mat');
%   stmdisp('lts2_levelling_vasteland.mat');
%
% Examples of plotting:
%
%   % cd to output folder
%   cd ../9_output
%
%   % Plot network map by epoch 
%   stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_ameland_v7_full_period.mat','saveplt','./');
%   stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_schier_v7_full_period.mat','saveplt','./');
%   stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_vasteland_v7_full_period.mat','saveplt','./');
%
%   % Estimate velocities and plot (using default 'refsystem'='min-velocities'
%   stmplotvel('../0_import_campaigns/lts2_levelling_ameland_v7_full_period.mat','saveplt','./');
%   stmplotvel('../0_import_campaigns/lts2_levelling_schier_v7_full_period.mat','saveplt','./');
%   stmplotvel('../0_import_campaigns/lts2_levelling_vasteland_v7_full_period.mat','saveplt','./');
%
%   % Plot the timeseries
%   stmplotseries('../0_import_campaigns/lts2_levelling_ameland_v7_full_period.mat','schranking','ssdnorm','marker',true,'saveplt','./');
%   stmplotseries('../0_import_campaigns/lts2_levelling_schier_v7_full_period.mat','schranking','ssdnorm','marker',true,'saveplt','./');
%   stmplotseries('../0_import_campaigns/lts2_levelling_vasteland_v7_full_period.mat','schranking','ssdnorm','marker',true,'maxPlots',4,'saveplt','./');
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
%   - you may ignore the close to singularity warnings in stmplotvel, if this bothers you, add the 
%     option (...,'ignoreStochModel',true)
%
%   - you may wish to explore various other schrankings options with stmplotseries
%
%        stmplotseries(...,'schranking','none', ...) no schranking at all, original spacetime matrix
%        stmplotseries(...,'schranking','ssdnorm',...) as in the main example, uses a two step procedure 
%                    minimizing differences between epochs and minimizing displacements for a stable area
%        stmplotseries(...,'schranking','ssdnorm','unstablearea','',...) as in the main example, but does
%                    not use a stable area 
%        stmplotseries(...,'schranking','ssdnorm','refpoints',{' ' , ... },'unstablearea','',...) as in 
%                   the main example, but uses reference points instead of polygon specifying a stable area 
%        stmplotseries(...,'schranking','ssdnorm','refpoints',{' ' , ... },...) augments the stable 
%                   area with reference points
%        stmplotseries(...,'schranking','min-velocities',...) use same method a stmplotvel (min-velocities in 
%                   stable area)
%
%      For instance to plot using a singe reference point as S-basis
%
%        stmplotseries('../0_import_campaigns/lts2_levelling_ameland.mat','schranking','ssdnorm','refPoints',{'000A2592'}, ...
%               'unstablearea','','marker',true);
%        stmplotseries('../0_import_campaigns/lts2_levelling_vasteland.mat','schranking','ssdnorm','refPoints','006B0021', ...
%               'unstablearea','','marker',true,'maxSubPlots',4,'maxPlots',6);
%
%      this ensures the displacements for the reference point will be zero, and for epochs without refenrence point the 
%      'schranking'='ssdnorm' method is used.
%
%      Also note that the 'saveplt' option is missing, so that the previous (default) plots will not be overwritten

% [End of script]



