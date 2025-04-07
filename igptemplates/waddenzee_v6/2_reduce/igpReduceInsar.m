%% Integrated geodetic processing - Reduce InSAR data
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project InSAR data reduction module.
%
% Input files
% - space time datasets with the InSAR data
% - file with points and epochs to reduce to
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
igpimport('mht');

%% Set the input filename and options

inputfilename = 'unreducedfiles.txt';

options = [];
options.verbose = 1;                 % Default verbosity level, higher is more output, 0 is almost nothing
options.doplots = 1;                 % Default plot level, 0 is no plots, higher is more detailed plotting

options.projectId = 'waddenzee';     % Project Id 
options.inputDir = '';               % Directory with the inputfiles
options.datasetId = '';              % DatasetId for the output data set, if empty, "_reduced" is added to the input datasetId.
options.outputTarget = 'create';     % Default mode for output STM {'create','overwrite','update'}, can be changed with flag

options.ROI = [];                    % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon (Default none)
options.POI = [-Inf +Inf];           % Period of interest [ dYearStart dYearEnd ] (Default none)

options.insarRadius = 500;            % [m], spatial search radius for spatial reduction
%options.diffEpochsMax = 0.2;          % [year], max abs epoch difference for temporal reduction [optional]
%options.diffEpochsMax = 0.1;          % [year], max abs epoch difference for temporal reduction [optional]
options.diffEpochsMax = [];          % [year], max abs epoch difference for temporal reduction [optional]
options.numEpochsMax = 5;           % maximum number of epochs for temporal reduction [optional]

options.splitDefoRegimes = true;    % Split in deformation regimes, true (default) or false

options.InSARModels = {... % InSAR models for selection of best model
            {'constant','linear'};...
            {'constant','linear','periodic'};...
           };
options.InSARModelsMinEpochs = [];   % Minimum number of epochs for partial model (breakpoint, heavyside)

options.stochModelParameterFile='../0_import_insar/insarStochasticModelParameters.txt'

%options.evaluationEpochs = './evaluation_epochs.mat';   % Name of the file with evaluation epochs (needed for reduction and epoch harmonization)
%options.evaluationPoints='./evaluation_points.mat';   % Name of the file with evaluation points (needed for harmonized point ids)


%% Call stmreduceinsar to do the reduction

inputfilenames = getinputfilenames(inputfilename,'insar'); % File with [insar] section
%inputfilenames = getinputfilenames(inputfilename); % File without [insar] section

%for k = 1:numel(inputfilenames)
%  outputfilename = [inputfilenames{k}(1:end-4) '_reduced.mat'];
%  stmreduceinsar(inputfilenames{k},outputfilename,options,'overwrite')
%end

for k=1:numel(inputfilenames)
   inputfilename=inputfilenames{k};
   [~,fname,fext]=fileparts(inputfilename);
   %outputfilename=strrep([fname fext],'decomposed','reduced');
   outputfilename=[fname '_reduced' fext];
   stmreduceinsar(inputfilename,outputfilename,options,'overwrite')
end

%% What to do next ...?
%
%  1. Inspect the stm using stmdisp or do plotting, for examples of 
%     plotting see./9_output/stmplot_examples.m 
%  2. Proceed with the integration

% [End of script]
