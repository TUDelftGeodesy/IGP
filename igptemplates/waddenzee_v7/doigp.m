%% Integrated geodetic processing - Waddenzee
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project levelling data import module.
%
% Main script

projdir=pwd;

%% Import campaign data 

cd 0_import_campaigns
run igpImportGPS
run igpImportLevelling
cd(projdir)

%% Import Insar Data

cd 0_import_insar
run igpImportInsar
cd(projdir)

%% Import GPS CORS data

cd 0_import_06gps
run igpImport06gps_202107
cd(projdir)

%% Decompose GPS CORS data

cd 1_decompose_06gps
run igpDecompose06gps
cd(projdir)

%% Reduce (select epochs and points)

cd 2_reduce
run igpSelect_waddenzee
cd(projdir)

%% Reduce (reduce datasets)

cd 2_reduce
run igpReduceCampaign
run igpReduceGnss        
run igpReduceInsar
cd(projdir)

%% Integrate

cd 3_integrate
run igpIntegrateWaddenzee2015
run igpIntegrateWaddenzee2009
run igpIntegrateWaddenzee2006
cd(projdir)

%% Predict

cd 4_predict
run igpPredictWaddenzee2015SimulatePointEpochList
run igpPredictWaddenzee2015
run igpExportWaddenzee2015

run igpPredictWaddenzeeLevellingVastelandSimulatePointEpochList
run igpPredictWaddenzeeLevellingVasteland
run igpExportWaddenzeeLevellingVasteland
cd(projdir)

