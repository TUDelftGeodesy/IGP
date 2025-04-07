%% Integrated geodetic processing - Reduce 06GPS GNSS data
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project GNSS data import module.
%
% Input files
% - space time datasets with the decomposed GNSS data from 06GPS
% - file with epochs to reduce to
%
% Ouputfile
% - reduced space time dataset
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

%% Set the input filenames, output filename and options

inputfilenames=getinputfilenames('unreducedfiles.txt','gnss');

roi=load('../2_reduce/roi_waddenzee.mat');

options=[];
options.verbose=1;                                  % Verbosity level
options.projectId='waddenzee';                      % Project Id that belongs to evaluation epochs and points (used to get the project file)
options.ROI=roi.sarROI;                             % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
options.rmafilt=21;                                 % Use 21 day moving average

%% Call stmreducegnss to do the reduction

% stmreducegnss(inputfilename,outputfilename,options,'overwrite')

for k=1:numel(inputfilenames)
   inputfilename=inputfilenames{k};
   [~,fname,fext]=fileparts(inputfilename);
   outputfilename=strrep([fname fext],'decomposed','reduced');
   stmreducegnss(inputfilename,outputfilename,options,'overwrite')
end

%% What to do next ...?
%
%  1. Inspect the stm using stmdisp or do plotting, for examples of 
%     plotting see./9_output/stmplot_examples.m 
%  2. Proceed with the integration

% [End of script]