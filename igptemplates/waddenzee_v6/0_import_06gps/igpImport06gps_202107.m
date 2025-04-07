%% Integrated geodetic processing - Import 06GPS GNSS data
%
% *Hans van der Marel*
%
% Integrated processing for NAM Integrated Geodetic Processing (IGP) 
% project GNSS data import module.
%
% Input files
% - GNSS data files from 06GPS
% - Meta data file with receiver and antenna names, steps, other events
% - Meteo data (downloaded from KNMI)
%
% Ouputfile
% - space time datasets with GNSS data
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

% Input file names and output file name

inputfiles='../../../igpdata/nam_gnss/06gps_nam_202107/*.csv';
outputfile='06gps_nam_202107.mat';

% Options for the GNSS conversion

options=[];
options.inputFormat='06gps';                               % Input format
options.stationAliases={ 'DZYL'    'DZY1'    ; ... 
                         'SCHI'    'SCH1'    ; ...
                         'HOOG'    'HOO9'    };            % Aliases
options.stationEvents='allEventsEdited.txt';               % Name of the file with receiver, antenna name, and steps

options.defaultAntennaType='TPSCR.G5        TPSH';         % IGRS stations have not been included in the station event file yet, so
options.defaultReceiverType='TPS NET-G5          ';        % we give them (temporarily) a default name

% Meteo options
%
% For the decomposition of the GPS signal we require temperature and air 
% pressure data. This data is available at the KNMI website
% <http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/uurgeg_280_2011-2020.zip>
%
% Name KNMI station     Number  Type    Latitude  Longitude
% ------------------    ------  ------  --------  ---------
% Eelde                 280     A/AWS   53 08     06 35 
% Nieuw Beerta          286     AWS     53 12     07 09      No pressure
% Leeuwarden            270     A/AWS   53 13     05 46 
% Hoorn (Terschelling)  251     AWS     53 23     05 21 
% Stavoren              267     AWS     52 53     05 23      No pressure
% Marknesse             273     AWS     52 42     05 53      No pressure 
% Lauwersoog            277     C/AWS   53 25     06 12      No pressure
% Hoogeveen             279     A/AWS   52 44     06 31 
%
% AWG-1                         PL/AWS  53 30     05 57      Cannot find data on KNMI site
%
% AWS=Automatic Weather Station, A=Airport, C=Coast, PL=Platform

% Structure with information on meteo stations

options.meteoProvider='KNMI';
% Name KNMI station     Number   Coordinates 
options.meteoStations = cell2struct({ ...
     'Eelde'                 280     [ 53+08/60  6+35/60 0 ] ; ... 
     'Nieuw Beerta'          286     [ 53+12/60  7+09/60 0 ] ; ...
     'Leeuwarden'            270     [ 53+13/60  5+46/60 0 ] ; ...
   %  'Hoorn (Terschelling)'  251     [ 53+23/60  5+21/60 0 ] ; ...
   %  'Stavoren'              267     [ 52+53/60  5+23/60 0 ] ; ...
   %  'Marknesse'             273     [ 52+42/60  5+53/60 0 ] ; ... 
   %  'Lauwersoog'            277     [ 53+25/60  6+12/60 0 ] ; ...
   %  'Hoogeveen'             279     [ 52+44/60  6+31/60 0 ] ; ...
   },{'name', 'number', 'coordinates' },2) ;
%options.meteoBaseURL ='http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/'; % Old http address not valid anymore, only https, but
%options.meteoBaseURL = 'https://www.knmi.nl/nederland-nu/klimatologie/uurgegevens/';         % https not supported by Matlab   
%options.meteoBaseURL = 'https://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/'; % Works end of 2024
options.meteoBaseURL = 'localhost';                                                           % Set to 'localhost' to skip automatic download, do this yourself
options.meteoFileTemplate = { 'uurgeg_###_2001-2010.txt' ; 'uurgeg_###_2011-2020.txt' ; 'uurgeg_###_2021-2030.txt' };

options.meteodir='../../../igpdata/nam_gnss/meteo_202107';    % Download meteo data into directory |./meteo|

%options.plotMeteo = true;          % Option to plot meteo data (default false)
%options.savePlots = true;          % Option to save plots as png (default false)

% Global attributes (use defaults defined by igpinit and modify)

globalAttrib=globalAttribDefault;   
globalAttrib.source = '06-GPS / Nederlandse Aardolie Maatschappij (NAM), The Netherlands.';
globalAttrib.technique = 'GNSS CORS';

options.globalAttrib=globalAttrib;

%% Do the conversion and save space time matrix to disk

[ts,st]=gnss2stm(inputfiles,outputfile,options,'overwrite');

%% What to do next ...?
%
%  1. Plot the network and/or timeseries
%  2. Decompose the timeseries using stmdecomposegnss 
%  2. Integrate the 2006+ levelling data, GPS campaign and InSAR 
%
% Examples of plotting:
%
%   cd ../9_output
%   stmplotmap('../0_import_06gps/06gps_nam_202107.mat')
%   stmplotmapbyepoch('../0_import_06gps/06gps_nam_202107.mat','saveplt','./')
%   stmplotseries('../0_import_06gps/06gps_nam_202107.mat','schranking','none','saveplt','./');
%
%   Notes: 
%
%   - because of the large number of epochs...
%
%        - stmplotmapbyepoch will not plot individual maps for each epoch
%        - stmplotvel will fail since stmvelocity cannot handle that many
%          epochs, anyway, 
% 
%               !!the decomposition will do a much better job!!
%
%     stmdecomposegnss also has a more elaborate trend model than lineair velocities
%
%   - plots are saved to pdf in the current folder, if you don't want this, leave out the 'saveplt' option, you can 
%     also specify another folder in 'saveplt'. The simple stmplotmap does not have a saveplt option
%
%   - if you want to run from another directory than ../9_ouput, add options to find the some of the files: 
%
%     stmplotmapbyepoch(...,'background','../9_output/igp-background.mat')
%     stmplotvel(...,'background','../9_output/igp-background.mat','unstablearea','../9_ouput/igp-unstablearea.mat')
%     stmplotseries(...,'unstablearea','../9_ouput/igp-unstablearea.mat');
%
%   - you may wish to explore some other options, see the build in Matlab help
%
% [End of script]
