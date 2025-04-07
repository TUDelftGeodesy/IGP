%% Integrated geodetic processing - Import InSAR data
%
% *Freek van Leijen, Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project InSAR data import module.
%
% Script to import InSAR datasets and create space-time matrices (one
% per dataset). Optionally, a classification of the measurement points
% and a separation in deformation regimes can be applied.
%
% Input files
% - InSAR data files
% - terrain height files [optional]
%
% Ouputfiles
% - space-time datasets with InSAR data
% - .txt file with a list of the generated space-time datasets


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

% Select input files for processing 

select='';               % Process all inputfiles 

% Input file names and definitions

insarInput = { ...
{ ...
'deformationFilename','lauwersmeer_s1_asc_gaussian90_deformation.csv', ...
'terrainFilename','lauwersmeer_s1_asc_gaussian90_ahn2.csv', ...
'headingAngle',350.0, ...
'system','s1', ...
'systemMode','IW', ...
'mode','Asc', ...
'polarization','VV' ...
'classificationOrder',{} ...
}, ...
{ ...
'deformationFilename','lauwersmeer_s1_dsc_gaussian90_deformation.csv', ...
'terrainFilename','lauwersmeer_s1_dsc_gaussian90_ahn2.csv', ...
'headingAngle',190.0, ...
'system','s1', ...
'systemMode','iw', ...
'mode','Dsc', ...
'polarization','VV', ...
'classificationOrder',{} ...
}, ...
{ ...
'deformationFilename','../nam_insar/radarsat2_desc_groningen_full_201812/U05_NAM_GTZH_U05_deformation.csv', ...
'terrainFilename','../nam_insar/radarsat2_desc_groningen_full_201812/U05_NAM_GTZH_U05_ahn2.csv', ...
'headingAngle',191.2, ...
'system','rsat2', ...
'systemMode','standard', ...
'mode','Dsc', ...
'polarization','' ...
}, ...
};

%% Set the options

% General options for the INSAR conversion

options=[];
options.inputFormat='Shell';                               % Input format
options.inputDir='../../../igpdata/nam_waddenzee/';        % Path to where InSAR files reside

% Global attributes (use defaults defined by igpinit and modify)

globalAttrib=globalAttribDefault;   
globalAttrib.source = 'SkyGeo / Nederlandse Aardolie Maatschappij (NAM), The Netherlands.';
globalAttrib.technique = 'INSAR';

options.globalAttrib=globalAttrib;

%% Do the conversion and save space time matrix to disk

datasets = [];

for k = 1:numel(insarInput)

   fprintf('\nProcessing InSAR file %d of %d ... \n',k,numel(insarInput));
   
   % Prepare inputfile name, outputfile name and options for insar2stm

   if strmatch('deformationFilename',insarInput{k}{1})
     inputfile = insarInput{k}{2};
   else
     error('The InSAR deformation file is not specified properly.');
   end
   [~,datasetId,~] = fileparts(inputfile);
   outputfile = [ datasetId '.mat' ];
   
   tmpoptions = options;
   for l=3:2:numel(insarInput{k})
      tmpoptions.(insarInput{k}{l}) = insarInput{k}{l+1};
   end

   % Call insar2stm to do the actual processing

   if contains(inputfile,select,'IgnoreCase',true')
      insar2stm({inputfile},outputfile,tmpoptions,'overwrite');
   else
      fprintf('Skip %s\n',inputfile)
   end
   
   % Save outputfile name

   datasets = [datasets cellstr(outputfile)]; 
end


%% Write file with output filenames  

fid = fopen('STMFilesInsar.txt','w');
numDatasets = numel(datasets);
for k = 1:numDatasets
  fprintf(fid,'%s\n',datasets{k});
end
fclose(fid);

%% What to do next ...?
%
%  1. Plot the network 
%  2. Integrate the 2006+ levelling data, GPS campaign and GPS CORS 
%
% Examples of plotting:
%
%   cd ../9_output
%   stmplotmap('../0_import_insar/U05_NAM_GTZH_U05_deformation.mat')
%   stmplotmap('../0_import_insar/lauwersmeer_s1_dsc_gaussian90_deformation.mat')
%   stmplotmap('../0_import_insar/lauwersmeer_s1_asc_gaussian90_deformation.mat')
%
%   Notes: 
%
%   - because of the large number of epochs...
%
%        - stmplotmapbyepoch will not plot individual maps for each epoch
%        - stmplotvel will fail since stmvelocity cannot handle that many
%          epochs
%
%   - because of the large number of points you have to be careful with
%     stmplotseries
%
%   - consider plotting after reducing the data in 2_reduce

% [End of script]
