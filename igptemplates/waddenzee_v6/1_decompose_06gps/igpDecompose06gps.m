%% Integrated geodetic processing - Decompose 06GPS GNSS data
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project GNSS data import module.
%
% Input file
% - space time datasets with tGNSS data from 06GPS (e.g created by igpImport06gps.m)
%
% Ouputfile
% - space time dataset with the decomposed time series
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
igpimport('tseries2');
igpimport('rdnaptrans');  % for plotting of maps

%% Set the input filenames, output filename and options

%inputfilename='06gps_nam_202101.mat';
%outputfilename='06gps_nam_decomposed_202101.mat';
inputfilename='../0_import_06gps/06gps_nam_202107.mat';
outputfilename='06gps_nam_202107_decomposed.mat';

options=[];
options.verbose=1;                                  % Verbosity level

options.doplots=0;                                  % Default plot level, 0 is no plots, higher is more detailed plotting

options.excludeFromCM={'AME2' 'AWG1' 'KOLH' 'NORG' 'NOR3' 'VEEN'};    % Stations to exclude from common mode and periodogram computations
options.includeHarmonics={ 'KOLH' 'NORG' 'NOR3' };           % Stations to include harmonic components in final estimate

options.removeCM=false;                             % If true, do a second iteration to remove the common mode effects (default false)

options.saveintermediate=true;                      % Save timeseries in tseries format

%options.globalAttrib=globalAttribDefault;                  % Use defaults defined by igpinit
options.globalAttrib=[];                                    % Copy global attributes from inputfile (is default)
options.globalAttrib.comment='Decomposed time series';      % Add comment attributes copied from NetCDF 


%% Call stmdecomposegnss to do the decomposition

stmdecomposegnss(inputfilename,outputfilename,options,'overwrite')

%% What to do next ...?
%
%  1. Plot the network and/or timeseries
%  2. Reduce the timeseries using stmselect and stmreducegnss, and integrate 
%     with levelling data, GPS campaign and InSAR 
%
% Examples of stmdisp and plotting:
%
%   stmdisp('../1_decompose_06gps/06gps_nam_202107_decomposed.mat')
%
%   cd ../9_output
%   stmplotmap('../1_decompose_06gps/06gps_nam_202107_decomposed.mat')
%   stmplotmapbyepoch('../1_decompose_06gps/06gps_nam_202107_decomposed.mat','saveplt','./')
%   stmplotseries('../1_decompose_06gps/06gps_nam_202107_decomposed.mat','schranking','none','saveplt','./');
%   stmplotseries('../1_decompose_06gps/06gps_nam_202107_decomposed.mat','schranking','none','item','auxData','saveplt','./')
%
%   Notes: 
%
%   - the results of the decomposition are save in the space time matrix,
%     see stmdisp output
%
%   - the first two plots are not any different than those from before the
%     decomposition
%
%   - because of the large number of epochs stmplotmapbyepoch will not plot 
%     individual maps for each epoch
%
%   - stmplotseries can be used to plot the resulting time series and to
%     plot the remainders (differences with the original) that is stored
%     in the aux data
%
%   - stmplotvel does not make sense for the decomposed result
%     (decomposition does a much better job)
%
%   - plots are saved to pdf in the current folder, if you don't want this, 
%     leave out the 'saveplt' option, you can also specify another folder 
%     in 'saveplt'. The simple stmplotmap does not have a saveplt option

% [End of script]