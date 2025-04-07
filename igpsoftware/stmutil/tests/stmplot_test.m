%% Space-time matrix (STM) plots  
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)

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

igpimport('stmutil');     % add all required toolboxes to the Matlab path
igpimport('crsutil');
igpimport('rdnaptrans');

%% Available space time matrices

projectdir = { ...
    '../' ; ...
    'd:/Surfdrive/Research/Projects/NAM/INTEGRATED PROCESSING/igpprojects/waddenzeehpc/' ...
             };

stmimport = { ...
   '0_import_campaigns/lts2_allgps_v6.0_20220111.mat' ; ...
   '0_import_campaigns/lts2_levelling_vasteland.mat' ; ...
   '0_import_campaigns/lts2_levelling_ameland.mat' ; ...
   '0_import_insar/lauwersmeer_s1_dsc_gaussian90_lauwersmeer_s1_dsc_gaussian90_deformation.mat' ; ...
   '0_import_insar/lauwersmeer_s1_asc_gaussian90_lauwersmeer_s1_asc_gaussian90_deformation.mat' };

stmproject = '2_reduce/waddenzee.mat';

stmreduced = { ...
   '2_reduce/06gps_nam_202107_reduced.mat' ; ...
   '2_reduce/lts2_allgps_v6.0_20220111_reduced.mat' ; ...
   '2_reduce/lts2_levelling_vasteland_reduced.mat' ; ...
   '2_reduce/lts2_levelling_ameland_reduced.mat' ; ...
   '2_reduce/lauwersmeer_s1_dsc_gaussian90_lauwersmeer_s1_dsc_gaussian90_deformation_reduced.mat' ; ...
   '2_reduce/lauwersmeer_s1_asc_gaussian90_lauwersmeer_s1_asc_gaussian90_deformation_reduced.mat' ; ...
   '2_reduce/U05_NAM_GTZH_U05_deformation_reduced.mat' };

stmintegrate = { ...
   '3_integrate/waddenzee_insarReducedWithDiag_LevDiag_test.mat' ; ...
   '3_integrate/waddenzee_insarReducedWithDiag_LevDiag.mat'; ...
   '3_integrate/waddenzee_insarReducedWithDiag_lp.mat' };


%% Select stmfile

%stmfile = fullfile(projectdir{1},stmimport{1})
%stmfile = fullfile(projectdir{1},stmproject{1})
stmfile = fullfile(projectdir{1},stmreduced{1})

%stmfile = fullfile(projectdir{1},stmintegrate{1})
%stmfile = fullfile(projectdir{2},stmintegrate{2})
%stmfile = fullfile(projectdir{2},stmintegrate{3})

%% Plot map with points

stmplotmap(stmfile)

%% Plot covariance matrix

stmplotcov(stmfile); %,'chol',true)

%% Plot map with evaluation points

stmfile1=fullfile(projectdir{2},stmproject);
stmfile2=fullfile(projectdir{2},stmintegrate{2});

stmplotprojectmap(stmfile1)
stmplotprojectmap(stmfile1,stmfile2)

%% Velocity plots (reduced datasets)

% Estimate velocities reduced dataset (single file - no printing)

stmplotvel(stmfile,'patch',true);

% Estimate velocities reduced dataset (single file - printing

stmplotvel(stmfile,'patch',true,'saveplt','./tmp');

%% Velocity plots (multiple datasets, with ROI)

roi=load('../2_reduce/roi_waddenzee.mat'); 
sarROI= roi.sarROI;
levROI= roi.levROI;

tmpfiles={};
for k=1:4
    tmpfiles{k}=fullfile(projectdir{1},stmreduced{k});
end
for k=5:numel(stmreduced)
    tmpfiles{k}=fullfile(projectdir{2},stmreduced{k});
end

stmplotvel(tmpfiles,'patch',true,'saveplt','./tmp','ROI',sarROI);

%% Velocity plots (integrated datasets) 

stmfile = fullfile(projectdir{1},stmintegrate{1})
%stmfile = fullfile(projectdir{2},stmintegrate{2})
%stmfile = fullfile(projectdir{2},stmintegrate{3})

stmplotvel(stmfile,'ignoreStochModel',true,'saveplt','./tmp');  % patch gives an error

