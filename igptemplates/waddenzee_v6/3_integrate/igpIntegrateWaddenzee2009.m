%% Integrated processing of geodetic observations - Waddenzee 2009 onwards
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)
% project using Waddenzee data.
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

%% Integrated processing 
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) project
%
% Syntax
% 
%    stmintegrate(inputfilenames,outputfilename,options)
%
% stmintegrate is called as script for debugging purposes (inputfilenames, 
% outputfile name and options must be defined in this script with the same
% name as used in stmintegrate
%
% Options can be specified as a structure 
%
%   options=[];
%   options.verbose=1;
%   options.pntCrdType='km/m';
%   options.aprioriDispl='simcase0c_truth.mat';
%   options.selectTechnique= { 'lev' 'gnss' };
%
% or cell array
%
%    options={'verbose',1,'pntCrdType','km/m'};
%
% It is also possible to select subsets of the techniques without redoing
% the simulation

%% Integrated processing for the Waddenzee - 2006.0 - 2020.5

inputfilenames='reducedfiles.txt';

options=[];
options.verbose=1;
options.plotCond=false;

options.selectTechnique= { 'insar' 'gnss' 'lev'};
options.filterLos=true;

roi=load('../2_reduce/roi_waddenzee.mat'); 
options.ROI= roi.sarROI;  
options.POI= [2009.3 Inf];

options.saveResiduals=true;

options.addStochModel.insar = { {'stdev(sigma=1)'} } ;                % add diagonal with 1mm2 
%options.addStochModel.lev = { {'houtenbos(s2t=0.16,pt=1.24)'} } ;    % add Houtenbos model with s2t=0.16 and pt=1.24
options.addStochModel.lev = { {'stdev(sigma=1)'} } ;                  % add diagonal with 1mm2 

% Force all GPS campaign data points to be included
%
% - for this option to work there must be for each epoch overlapping points with other datasets
% - this measn, the last epoch, which does not have the CORS data, must be 
%   removed; hence the redefined POI

options.forceDatasetPoints='lts2_allgps';     
options.POI= [2009.3 2020.5];

outputfilename='waddenzee_2009.mat';

stmintegrate(inputfilenames,outputfilename,options,'overwrite')

%% What to do next ...?
%
%  1. Inspect the stm using stmdisp
%  2. Plot the residuals 
%  3. Display the time series
%  4. Proceed with the prediction
%
% Plot residuals (interactive display)
%
%    stmresplot(outputfilename)
%
% For examples of other plotting see./9_output/stmplot_examples.m 
  
% [End of script]
