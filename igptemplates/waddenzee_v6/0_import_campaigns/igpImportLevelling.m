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

inputfile={'../../../igpdata/nam_waddenzee/lts2_alllevelling_gpsbaseline_v5.2_20210323.nc'}; %format {'file1.nc';file2.nc';...}

% OVERVIEW OF OPTIONS FOR SELECTING NETWORKS
% ------------------------------------------

% All campaigns

% options.includePrjName = { };
% options.excludePrjName={ };
% outputfile='lts2_levelling_all.mat';

% Hydrostatisch ?
%
% 279H05                1986.476  1986-Jun-24   1819     48
% 344W01                1993.975  1993-Dec-23    178     31
% 371W00                1998.329  1998-May-01     91    195
% 373W31                2001.504  2001-Jul-04    299     59
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

% Ameland
%
% 289W05                1987.837  1987-Nov-03     94     73
% 289W16                1988.749  1988-Oct-01    333     45
% 289W20                1990.085  1990-Feb-01    488     48
% 289W26                1991.084  1991-Feb-01    262     46
% 289W34                1992.045  1992-Jan-18    249     48
% 279W22                1992.371  1992-May-16    119      5
% 332W02                1993.086  1993-Feb-01    261     74
% 289W37                1994.088  1994-Feb-02     41     45
% 342W04                1995.060  1995-Jan-23    355     45
% 342W05                1996.029  1996-Jan-12    354     43
% 342W10                1997.025  1997-Jan-10    257     42
% options.includePrjName = { ... 
%   '289W05'  '289W16'  '289W20'  '289W26'  '289W34'  '279W22'  '289W37'  '332W02'  '342W04'  '342W05'  '342W10' ... 
%   'aml1998'  'aml1999' 'aml2003'  'AMEL0409' 'aml2005' ...
%  '70604-001' '231-70604' '231-70812' 'Amel_2009' 'NAM_AM2011' 'NAM_AM2014' 'NAM_AM2017' 'NAM_AM2020' ...
% };
% options.excludePrjName={ };
% outputfile='lts2_levelling_ameland.mat';

% Vasteland
%
% 231                   1981.496  1981-Jul-01      0    135
% 300                   1987.580  1987-Aug-01    403    475
% 289W22                1990.367  1990-May-15    103     75
% 289W31                1991.363  1991-May-14    102     74
% 332W04                1993.488  1993-Jun-28     13    471
% 342W08                1996.322  1996-Apr-28    107    300
% 342W12                1997.447  1997-Jun-13    154    360
% 364W00                1998.424  1998-Jun-05     35    579
% 365W37-SEC            1999.749  1999-Oct-02    234    302 
% 370W27                2000.686  2000-Sep-08     24    169
% GRON0306              2003.457  2003-Jun-17    135    634
% options.includePrjName = { ...
%   '231'  '300'  '289W22'  '289W31'  '332W04'  '342W08'  '342W12'  '364W00'  '365W37-SEC'  '370W27'  'GRON0306' ...
%   'NAM_GF2008' 'NAM_LM2011' 'NAM_LM2015'  'NAM_MA2006' 'RWS_GF2013'  'RWS_GF2018' ...
% };
% options.excludePrjName={ '231-'};
% outputfile='lts2_levelling_vasteland.mat';

% Options for ROI
% options.POI=[2006,+Inf];
% options.includePrjName={ '231-70604'  '231-70812'  '70604-001' 'Amel_2009' 'NAM_AM2011' 'NAM_AM2014' 'NAM_AM2017' 'NAM_AM2020'};
% options.includePrjName={'NAM_GF2008' 'NAM_LM2011' 'NAM_LM2015'  'NAM_MA2006' 'RWS_GF2013'  'RWS_GF2018' };


%% All - full period

outputfile='lts2_levelling_v6_all.mat';

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

outputfile='lts2_levelling_vasteland_full_period.mat';

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

outputfile='lts2_levelling_ameland_full_period.mat';

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
%   stmdisp('lts2_levelling_vasteland.mat');
%
% Examples of plotting:
%
%   % cd to output folder
%   cd ../9_output
%
%   % Plot network map by epoch 
%   stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_ameland.mat','saveplt','./');
%   stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_vasteland.mat','saveplt','./');
%
%   % Estimate velocities and plot (using default 'refsystem'='min-velocities'
%   stmplotvel('../0_import_campaigns/lts2_levelling_ameland.mat','saveplt','./');
%   stmplotvel('../0_import_campaigns/lts2_levelling_vasteland.mat','saveplt','./');
%
%   % Plot the timeseries
%   stmplotseries('../0_import_campaigns/lts2_levelling_ameland.mat','schranking','ssdnorm','marker',true,'saveplt','./');
%   stmplotseries('../0_import_campaigns/lts2_levelling_vasteland.mat','schranking','ssdnorm','marker',true,'maxPlots',4,'saveplt','./');
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



