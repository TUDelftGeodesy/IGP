%% Integrated geodetic processing - Reduce datasets
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project evaluation point and epoch selection.
%
% Input files
% - file with list of original space time data files
% - original (non reduced) space time data file
%
% Ouputfile
% - projectFile with selected evaluation points and epochs
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

%% Set the input filenames, output filename and options

% Input file names and output file name

inputfilenames = 'unreducedfiles.txt';
outputfilename = 'waddenzee.mat';    % Name of the projectFile

% Options for the evaluation point and epoch selection
%
% Commonly used options are
%
%     projectId=''                Project Id (if empty, then the outputfilename is used as projectId) 
%
%     selectTechnique= { ... }    Cell array with techniques to include (default is {'lev' 'gnss' 'insar'})
%     selectDataset= ''           Cell array with datasetIds to select, you can give a pattern to select, default is all
%
%     selectEpochs={ ... }        Cell array with techniques or modes used for selecting epochs (remaining techniques 
%                                 will be used for desification), default is {'lev' 'grav' 'Campaign') 
%     temporalTolerance = 0.1     Time tolerance [decimal years] for selection
%     numEpochsYear = 1;          Number of epochs per year for densification
%
%     selectPoints={ ... }        Cell array with techniques or modes used for selecting points (remaining techniques
%                                 will be used for densification), default is {'lev' 'gnss' 'grav'}
%     spatialTolerance = 0.1      Spatial tolerance [km] for merging points in selection
%     numPointsKm2 = 0.5          Density for densification [km2]
%     spatialMethod = 'grid'      Densification method, only 'grid' for now
%
%     globalAttributes=[]         Structure with global attributes 

options = [];
options.temporalTolerance = 1/12;   % time tolerance [decimal years] for selection
%options.numEpochsYear = 1;        % number of epochs per year for densification

%options.spatialTolerance = 0.1;    % spatial tolerance [km] for merging points in selection
options.spatialTolerance = 0.5;    % spatial tolerance [km] for merging points in selection
options.numPointsKm2 = 0.5;        % density for densification [km2]
options.spatialMethod = 'grid';    % densification method, only 'grid' for now

roi=load('../2_reduce/roi_waddenzee.mat');
sarROI =roi.sarROI;

options.densificationROI=sarROI;
%options.densificationROI={'lev'};
%options.densificationROI={'lev','gnss'};
%options.densificationROI='U22';
%options.shrink=.8;
%options.densificationROI= [53.25 6.5 ; 53.46 6.9 ];

options.verbose=0;
options.doplots=0;

options.globalAttrib=globalAttribDefault;   % Global Attributes from igpinit 

%% Select evaluation epochs and points, and write projectFile 

stmselect(inputfilenames,outputfilename,options,'overwrite');

%% What to do next ...?
%
%  1. Plot the network with evaluation points
%  2. Reduce the datasets and proceed with the integration
%
% Examples of plotting:
%
%   cd ../9_output
%   stmplotprojectmap('../2_reduce/waddenzee.mat')
%   stmplotprojectmap('../2_reduce/waddenzee.mat','../3_integrate/waddenzee_2006.mat')
%   exportgraphics(gca, 'evaluation_points_2006.pdf', 'ContentType', 'vector')
%
%   Notes: 
%
%   - the second example shows only the points that are actually used in
%     the integration
%
%   - there is not an option to save plots, you have to do this yourself 
%
%   - see also ./9_output/stmplot_examples.m 

% [End of script]
