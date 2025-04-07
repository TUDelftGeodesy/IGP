%% Space-time matrix (STM) plots  
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP)
%
% ----------------------------------------------------------------------------------
% Stmutil contains several top-level functions printing and plotting of space 
% time matrices:
% 
%    stmdisp           - Display basic information on space time matrix datasets 
%    stmdiff           - Diff (compare) two space time matrix datasets
%
%    stmplotmap        - Plot map with the points in the stm dataset 
%    stmplotprojectmap - Plot map of the evaluation points 
%    stmplotseries     - Plot timeseries from space time matrix datasets 
%    stmplotcov        - Plot space time matrix covariance matrix
%    stmplotvel        - Plot velocities estimated by stmvelocity 
%    stmresplot        - Plot residuals and test statistics from stmintegrate
%
% The top-level functions operate on space time matrices directly and are typically
% used to print and plot information from space time matrix datasets from the
% command line or scripts. 
% ----------------------------------------------------------------------------------
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

igpimport('stmutil');     % add all required toolboxes to the Matlab pathopt.sreg
igpimport('crsutil');
igpimport('rdnaptrans');

cd ../9_output

%% GPS CORS data

% plot network map (stmplotmapbyepoch will not plot individual maps for each epoch)
stmplotmap('../0_import_06gps/06gps_nam_202312.mat')
stmplotmapbyepoch('../0_import_06gps/06gps_nam_202312.mat','saveplt','./plt_imports/')
% Plot the timeseries
stmplotseries('../0_import_06gps/06gps_nam_202312.mat','schranking','none','saveplt','./plt_imports/');

% results from the decomposition
stmplotseries('../1_decompose_06gps/06gps_nam_202312_decomposed.mat','schranking','none','saveplt','./plt_imports/');
stmplotseries('../1_decompose_06gps/06gps_nam_202312_decomposed.mat','schranking','none','item','auxData','saveplt','./plt_imports/')


%% GPS Campaign data

% Plot network map by epoch -> plots are saved to a pdf 
stmplotmapbyepoch('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','saveplt','./plt_imports/');
% Estimate velocities and plot (using default 'refsystem'='min-velocities')
stmplotvel('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','saveplt','./plt_imports/');
% Plot the timeseries
stmplotseries('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','schranking','none','marker',true,'saveplt','./plt_imports/');
stmplotseries('../0_import_campaigns/lts2_allgps_v7.0_20230116.mat','schranking','ssdnorm','marker',true);

%% Levelling data (full period)

% Plot network map by epoch
stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_v7_all.mat','maxepochs',60,'saveplt','./plt_imports/');
stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_ameland_v7_full_period.mat','saveplt','./plt_imports/');
stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_schier_v7_full_period.mat','saveplt','./plt_imports/');
stmplotmapbyepoch('../0_import_campaigns/lts2_levelling_vasteland_v7_full_period.mat','saveplt','./plt_imports/');

% Estimate velocities and plot (using default 'refsystem'='min-velocities')
stmplotvel('../0_import_campaigns/lts2_levelling_ameland_v7_full_period.mat','saveplt','./plt_imports/');
stmplotvel('../0_import_campaigns/lts2_levelling_schier_v7_full_period.mat','saveplt','./plt_imports/');
stmplotvel('../0_import_campaigns/lts2_levelling_vasteland_v7_full_period.mat','saveplt','./plt_imports/');

% Plot the timeseries
stmplotseries('../0_import_campaigns/lts2_levelling_ameland_v7_full_period.mat','schranking','ssdnorm','marker',true,'saveplt','./plt_imports/');
stmplotseries('../0_import_campaigns/lts2_levelling_schier_v7_full_period.mat','schranking','ssdnorm','marker',true,'saveplt','./plt_imports/');
stmplotseries('../0_import_campaigns/lts2_levelling_vasteland_v7_full_period.mat','schranking','ssdnorm','marker',true,'maxPlots',4,'saveplt','./plt_imports/');

%% InSAR data (only network maps)

stmplotmap('../0_import_insar/U05_NAM_GTZH_U05_deformation.mat')
exportgraphics(gca,'./plt_imports/U05_NAM_GTZH_U05_deformation.pdf','Resolution',300);
stmplotmap('../0_import_insar/lauwersmeer_s1_dsc_gaussian90_deformation.mat')
exportgraphics(gca,'./plt_imports/lauwersmeer_s1_dsc_gaussian90_deformation.pdf','Resolution',300);
stmplotmap('../0_import_insar/lauwersmeer_s1_asc_gaussian90_deformation.mat')
exportgraphics(gca,'./plt_imports/lauwersmeer_s1_asc_gaussian90_deformation.pdf','Resolution',300);

%% Plot map with evaluation points (for reduced datasets and integration)

stmproject='../2_reduce/waddenzee.mat';
stmintegrate =  { ...
     '../3_integrate/waddenzee_2006.mat' ; ...
     '../3_integrate/waddenzee_2009.mat' ; ...
     '../3_integrate/waddenzee_2015.mat' };

stmplotprojectmap(stmproject)

stmplotprojectmap(stmproject,'../3_integrate/waddenzee_2006.mat')
exportgraphics(gca, './plt_reduced/evaluation_points_2006.pdf', 'ContentType', 'vector')

stmplotprojectmap(stmproject,'../3_integrate/waddenzee_2009.mat')
exportgraphics(gca, './plt_reduced/evaluations_point_2009.pdf', 'ContentType', 'vector')

stmplotprojectmap(stmproject,'../3_integrate/waddenzee_2015.mat')
exportgraphics(gca, './plt_reduced/evaluations_point_2015.pdf', 'ContentType', 'vector')

%% Velocity plots and time series plots for reduced datasets 

% get the paths to all reduced datasets
stmreduced = getinputfilenames('../3_integrate/reducedfiles.txt');
roi=load('../2_reduce/roi_waddenzee.mat'); 
sarROI= roi.sarROI;

% plot the velocities for all datasets and create three pdf with appended results (no loop needed)
stmplotvel(stmreduced,'ROI',sarROI,'saveplt','./plt_reduced/','append',true)

%% plot the time series for reduced datasets

for k=1:numel(stmreduced)
    if contains(stmreduced{k},'_s1_') 
        % the S1 dataset still contain "empty" shallowlos layers, we don't want to plot them
        stmplotseries(stmreduced{k},'schranking','ssdnorm','marker',true,'maxPlots',10,'types','los','saveplt','./plt_reduced/')
    else
        stmplotseries(stmreduced{k},'schranking','ssdnorm','marker',true,'maxPlots',10,'saveplt','./plt_reduced/')
    end
end


%% Residual plots for integrated datasets
%
%  Simple statement:  
% 
%      stmresplot(integrated_datasetname)
%
%  Complicating factors
%      we have multiple datasets
%      stmresplot has not yet an option to save plots, so we have to do
%      this ourselves

stmintegrate =  { ...
     '../3_integrate/waddenzee_2015.mat' ; ...
     '../3_integrate/waddenzee_2009.mat' ; ...
     '../3_integrate/waddenzee_2015.mat' };

for l=1:numel(stmintegrate)
   % stmresplot cannot save plots, so we do this outside
   hplt=findall(groot,'Type','figure');
   numfig=numel(hplt);
   [~,pdfname]=fileparts(stmintegrate{l});
   stmresplot(stmintegrate{l})
   hplt=findall(groot,'Type','figure');
   for k=numfig+1:numel(hplt)
       exportgraphics(hplt(k), ['./plt_integrate/' pdfname '_residuals.pdf'], 'ContentType', 'vector','Append',true)
   end
end

%% Velocity plots for integrated datasets 

% 2006 dataset is too big for my laptop ... try 2015 
% for l=1:numel(stmintegrate)
%    %stmplotvel(stmintegrate{l},'ignoreStochModel','true','saveplt','./plt_integrate');
%    stmplotvel(stmintegrate{l},'saveplt','./plt_integrate');
% end
stmplotvel('../3_integrate/waddenzee_2015.mat','saveplt','./plt_integrate');

%% Time series plots for integrated datasets

for l=1:numel(stmintegrate)
   stmplotseries(stmintegrate{l},'schranking','ssdnorm','marker',true,'maxPlots',30,'saveplt','./plt_integrate/')
   % resulting pdf is always 'integrated_obsData.pdf, give it a more meaningful name
   [~,pdfname]=fileparts(stmintegrate{l});
   movefile('./plt_integrate/integrated_obsData.pdf',['./plt_integrate/' pdfname '_series.pdf'])
end


