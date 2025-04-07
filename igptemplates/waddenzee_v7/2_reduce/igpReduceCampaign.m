%% Integrated geodetic processing - Reduce campaign data
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project campaign reduce module
%
% Input files
% - space time dataset with the campaign data to reduce
% - project file
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

%% Set the input filenames, output filename and options

% inputfilename='../0_import_levelling/gron_levelling_flaggedOutliers_v2.3.02_20191030.mat';
% outputfilename='gron_levelling_flaggedOutliers_v2.3.02_20191030_reduced.mat';

inputfilenames=getinputfilenames('unreducedfiles.txt','campaign');

options=[];
options.verbose=1;                                  % Verbosity level
options.projectId='waddenzee';                      % Project Id that belongs to evaluation epochs and points

%% Call stmreducecampaign to do the reduction

% stmreducecampaign(inputfilename,outputfilename,options,'overwrite')

for k=1:numel(inputfilenames)
   inputfilename=inputfilenames{k};
   [~,fname,fext]=fileparts(inputfilename);
   outputfilename=[ fname '_reduced' fext];
   stmreducecampaign(inputfilename,outputfilename,options,'overwrite')
end

%% What to do next ...?
%
%  1. Inspect the stm using stmdisp or do plotting, for examples of 
%     plotting see./9_output/stmplot_examples.m 
%  2. Proceed with the integration

% [End of script]