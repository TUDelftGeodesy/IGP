function [ts,st]=gnss2stm(inputfilenames,outputfilename,varargin)
%gnss2stm   Import gnss data files into space time matrix dataset
%   GNSS2STM(INPUTFILENAMES,OUTPUTFILENAME,OPTIONS) imports GNSS data files
%   given by INPUTFILENAMES into a single Space Time Matrix dataset with 
%   the name OUTPUTFILENAME. INPUTFILENAMES must be the name of an file 
%   containing the input filenames, a character string with wildcards, 
%   or a cell array with the input filenames. OUTPUTFILENAME is a character 
%   string with the name of the output Space Time Matrix dataset. OPTIONS 
%   is a cell array or structure with the processing options, if empty, 
%   or when fields are missing, default values are used.
%
%   GNSS2STM(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   STAT=GNSS2STM(...) returns a status code STAT. A status code of 0
%   indicates success.
%
%   OPTIONS is a cell array or structure with the processing options, 
%   if empty, or when fields are missing, default values are used. Valid 
%   OPTIONS are
%
%      verbose=0                Verbosity level, higher is more output, 0 is almost nothing
%
%      inputDir=''              Directory with the inputfiles
%      inputFormat='06gps'      Input format (Default '06gps' csv files)
%
%      stationAliases={}        n-by-2 cell array {'old' new' ; ...} with station aliases (Default none)
%      stationEvents=''         Name of the file with the station receiver name, antenna name, steps and other events
%      defaultAntennaType='UNKNOWN             '     
%                               Default antenna type in case station is not found in the events file
%      defaultReceiverType='UNKNOWN             '    
%                               Default receiver type in case station is not found in the events file
%
%      meteoProvider=''         Meteo data provider, default none, currently only KNMI supported 
%      meteoStations=struct()   Structure with meteo station data (format is provider dependent)
%      meteoBaseURL=''          Base URL where to find the data (or 'localhost')
%      meteoFileTemplate={}     Character cell array with meteo data filename templates
%      meteodir='meteo'         Directory where to store the meteo files (Default meteo)
%
%      plotMeteo=false          Option to plot meteo data (Default false)
%      savePlots=false          Option to save plots to disk (Default false)
%
%      ROI=[]                   Region of interest, as [latmin lonmin ; latmax lonmax] bounding 
%                               box, or lat/lon polygon, or kml/shape file (Default all) *)
%      POI=[-Inf +Inf]          Period of interest [ dYearStart dYearEnd ] or 
%                               [ dYearStart dYearEnd ; dYearStart dYearEnd ;...  ]  (Default all) *)
%
%      globalAttrib=[]          Struct with global attributes for output dataset (Empty by default) 
%
%      projectId=''             Project Id (default directory of the outputfile, or empty if no directory is specified)
%
%   [Note: ROI and POI are NOT implemented at the moment...]
% 
%   Examples:
%      gnss2stm('csv/*.csv','groningen_GNSS.mat')      
%      gnss2stm('csv/*.csv','groningen_GNSS.mat',options,'update')      
%
%   See also STM, STMCHECKARGUMENTS, STMREAD and STMDIFF.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:   2 November 2020 by Hans van der Marel
% Modified: 11 March 2021 by Hans van der Marel
%              - added option to read meteo data and save this as epoch
%                attributes, supports KNMI meteo data download and reading
%              - added option to set default receiver and antenna name
%           29 August 2021 by Hans van der Marel
%              - Changes related to new URLs for the meteo data, if url
%                is 'localhost' no meteo data is downloaded, but assumed to
%                be present in the meteo data directory already. This 
%                became necessary since Matlab does not support https and
%                KMNI does not support http anymore. Manual downloads are
%                now required.
%           22 September 2021 by Hans van der Marel
%              - Save options as datasetAttributes and meteostation names as
%                techniqueAttributes in output space time matrix
%              - Added options to disable plotting meteo data (default no plots)
%           24 October 2021 by Hans van der Marel
%              - removed local function bbox2poly
%              - commented out ROI and POI options (not yet implemented)
%              - added global attributes to output stm
%              - added description of options to the help

% TO DO:
%              - implement ROI and POI (use getPointMask and getEpochMask)


%% Check the input arguments and options

if nargin < 2
   error('This function expects at least two input arguments.')
end

progname='gnss2stm';

% Default options

opt.verbose=0;                                  % Default verbosity level, higher is more output, 0 is almost nothing

opt.inputDir='';                                % Directory with the inputfiles
opt.inputFormat='06gps';                        % Input format

opt.stationAliases={};                          % Station aliases
opt.stationEvents='';                           % Name of the file with the station receiver name, antenna name, steps and other events
opt.defaultAntennaType='UNKNOWN             ';  % Default antenna type in case station is not found in the events file
opt.defaultReceiverType='UNKNOWN             '; % Default receiver type in case station is not found in the events file

opt.meteoProvider='';                           % Meteo data provider, default none, currently only KNMI supported 
opt.meteoStations=struct([]);                   % Structure with meteo station data (format is provider dependent)
opt.meteoBaseURL='';                            % Base URL where to find the data (or 'localhost')
opt.meteoFileTemplate={};                       % Character cell array with meteo data filename templates
opt.meteodir='meteo';                           % Directory where to store the meteo files (Default meteo)

opt.plotMeteo=false;                            % Option to plot meteo data (Default false)
opt.savePlots=false;                            % Option to save plots to disk (Default false)

%opt.ROI=[];                                     % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon, kml/shape file (Default none)
%opt.POI=[-Inf +Inf];                            % Period of interest [ dYearStart dYearEnd ] or [ dYearStart dYearEnd ; dYearStart dYearEnd ; ... ] (Default none)

opt.globalAttrib=struct([]);                     % Struct, if empty (default) globalAttributes in stm is empty

opt.projectId='';                               % Project Id (default directory of the outputfile, or empty if no directory is specified)

% Duplicate output to file, and catch errors, and start timing

try

[~,outputfileroot]=fileparts(outputfilename);
diary([ outputfileroot '_' datestr(now,30) '.log' ])

fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Handle wildcards

% if ischar(inputfilenames)
%    d=dir(inputfilenames);
%    if numel(d) > 1
%      inputfilenames={d.name};
%      opt.inputDir=unique({d.folder});
%    end
% end

% Check the options and if necessary overwrite the default values

[inputfilenames,outputfilename,opt]= ...
    stmcheckarguments(inputfilenames,outputfilename,opt,varargin{:});
if isempty(outputfilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    return;
end


%% Prepare for the import

% Set constants

[a,f] = inqell('WGS-84');
e2 = 2*f - f^2;

oneday=1/365.25;

% Read the event file

% eventfile='d:\Surfdrive\Research\Projects\NAM\GNSS_APM\datasets\allEventsEdited.txt';
% for k=1:numel(ts)
%  ts(k).events=getEvents(eventfile,ts(k).station,ts(k).year(1),Inf);
% end

if ~isempty(opt.stationEvents)
   if exist(opt.stationEvents,'file')
      allEvents=readtable(opt.stationEvents,'format','%s %s %s %s %s');
   else
      error('Incorrect stationEvent file, please specify the correct name.')     
   end
else
   warning('Missing stationEvent file, receiver and antenna changes, and other events, are ignored. Proceed with caution!!')
   allEvents=table();
end

% This is how to use the allEvent table 
%
% eventtable=table2struct(allEvents(ismember(allEvents.station,station),2:end));
% for k=1:size(eventtable,1)
%   events(k).year=date2dyear(eventtable(k).date,'yyyy-mm-ddTHH:MMZ');
%   events(k).type=eventtable(k).type;
%   events(k).name=eventtable(k).name;
% end
% events=trimEvents(events,tfirst,tlast,'station',station);

% Set aliases (example)
%
%               06-GPS    new name
% opt.aliases={ 'DZYL'    'DZY1'    ; ...
%               'SCHI'    'SCH1'    ; ...
%               'HOOG'    'HOO9'    };

aliases=opt.stationAliases;


%% Read the input files into a time series structure

switch opt.inputFormat

    case '06gps'

        % Read 06-GPS comma separated files 
        %
        % Missing and duplicated data is removed from the csv files using the 
        % function namnodata.m
        %
        % The hourly samples in the csv files are downsampled to daily samples
        % using the mean over the interval.

        softwareName='Geo++';
        
        ts=[];
        for k=1:numel(inputfilenames)

           %staid=strrep(staid,'.csv','');
           [~,staid]=fileparts(inputfilenames{k});
           fprintf('\nImporting %s ... \n\n',staid)

           % Read csv file

           fid=fopen(fullfile(opt.inputDir,inputfilenames{k}));
           %2006-06-02 00:00,53,22,15.04169,6,9,8.59169,45.2852
           c=textscan(fid,'%s %d %d %f %d %d %f %f','delimiter',',');
           fclose(fid);

           epoch=datenum(char(c{1}),'yyyy-mm-dd HH:MM');
           lat=double(c{2})+double(c{3})/60+c{4}/3600;
           lon=double(c{5})+double(c{6})/60+c{7}/3600;
           hgt=c{8};
           plh=[lat*pi/180 lon*pi/180 hgt];

           plh0=mean(plh);
           N = a ./ sqrt(1 - e2 .* sin(plh0(1)).^2);
           M = N * (1 -e2) / ( 1 - e2 .* sin(plh0(1)).^2 );
           flat= M;
           flon= N*cos(plh0(1));

           neu=[ (plh(:,1)-plh0(1))*flat   (plh(:,2)-plh0(2))*flon   plh(:,3)-plh0(3) ];
           year=date2dyear(epoch);

           % Remove periods without data

           nodata=namnodata(epoch,neu,0.0001);

           plh(nodata,:)=nan;
           neu(nodata,:)=nan;

           % Detect and remove duplicate epochs

           duplicates=find(diff(epoch) <=0 );
           if ~isempty(duplicates)
             fprintf('The data contains duplicate or decreasing epochs, remove...\n')
             for k1=1:numel(duplicates)
               kk=duplicates(k1)+1;
               fprintf('%s  %s   %10.4f %10.4f %10.4f    %10.4f %10.4f %10.4f\n',datestr(epoch(kk-1),31),datestr(epoch(kk),31),neu(kk-1,:),neu(kk,:))
             end
             epoch(duplicates)=[];
             year(duplicates)=[];
             plh(duplicates,:)=[];
             neu(duplicates,:)=[];
           end

           % Downsample to daily estimates (using mean or median)

           xq=floor(min(epoch))+.5:1:ceil(max(epoch))-0.5;
           %[yqmm,syqmm,nyqmm]=mminterp1(epoch,neu,xq);
           [yqma,syqma,nyqma]=mainterp1(epoch,neu,xq);

           idx=nyqma >= 20 & ~any(isnan(yqma),2);

           epoch=xq(idx);
           epoch=epoch(:);
           neu=yqma(idx,:);
           sneu=syqma(idx,:);

           year=date2dyear(epoch); 

           % save data in structure

           ts(k).station=upper(staid);     % station name (4 letter abbreviation, upper case)
           ts(k).plh0=plh0;                % ref. latitude, longitude and height ([deg] and [m])
           ts(k).vneu0=zeros(1,3);         % a-priori velocity [m/y]
           ts(k).t0=round(mean(year));     % reference time in years (must be rounded to whole years!!)

           ts(k).events=struct([]);        % structure array with events will be added later

           ts(k).epoch=epoch;              % epoch as Matlab datenumber
           ts(k).year=year;                % epoch in decimal years 
           ts(k).neu=neu;                  % neu time series
           ts(k).sneu=0.001*ones(size(neu)); % a-priori standard deviation

           ts(k).geo2m=[flat flon 1];      % conversion factors for geographic coordinates [deg] to meters
           ts(k).plh=plh;

           ts(k).filename=fullfile(opt.inputDir,inputfilenames{k});
           
        end
        
    otherwise
        
        error(['Unsupported data format ' opt.inputFormat ])        

end

%% Rename stations and add event data to time series structure

numAliases=0;
for k=1:numel(ts)
    % check aliases for station name and make it uppercase
    staid=ts(k).station;
    if any(ismember(aliases(:,1),staid))
       staidold=staid;
       staid=char(aliases(ismember(aliases(:,1),staid),2));
       fprintf('Station %s has been renamed, new name is %s\n',staidold,staid)
       ts(k).station=staid;
       numAliases=numAliases+1;
    end
end
if numAliases ~=0
   fprintf('\n%d stations have been renamed\n\n',numAliases)
end

numMissingEvents=0;
for k=1:numel(ts)
    % get event data from the global event table
    staid=ts(k).station;
    year=ts(k).year;
    eventtable=table2struct(allEvents(ismember(allEvents.station,staid),2:end));
    events=struct([]);
    if numel(eventtable) > 0
       for l=1:size(eventtable,1)
          events(l).year=date2dyear(eventtable(l).date,'yyyy-mm-ddTHH:MMZ');
          events(l).type=eventtable(l).type;
          events(l).name=eventtable(l).name;
       end
       events=trimEvents(events,min(year)-oneday/2,max(year)+oneday/2,'station',staid);
    else
       fprintf('Station %s has no meta data (events) defined (receiver type set to %s, antenna type set to %s).\n',staid,opt.defaultReceiverType,opt.defaultAntennaType)
       events(1).year=min(year)-oneday/2;
       events(1).type='ANT';
       events(1).name=opt.defaultAntennaType;
       events(2).year=min(year)-oneday/2;
       events(2).type='REC';
       events(2).name=opt.defaultReceiverType;
       numMissingEvents=numMissingEvents+1;
    end
    ts(k).events=events;
end
if numMissingEvents ~=0
   fprintf('\n%d stations have no receiver or antenna type data, proceed with caution!!\n\n',numMissingEvents)
end

pntEvents=cell(numel(ts),1);
pntDyearRange=nan(numel(ts),2);
for k=1:numel(ts)
   events=ts(k).events;
   c={};
   for l=1:numel(events)
       c{l}=strjoin({ events(l).type, datestr(dyear2date(events(l).year),30), events(l).name },',');
   end
   pntEvents{k}=strjoin(c,';');
   pntDyearRange(k,:)=[  min(ts(k).year)-oneday/2  max(ts(k).year)+oneday/2 ];
end

%% Convert time series structure into space time matrix format

numPoints=numel(ts);

sites={ts(:).station};
plh0=[ts(:).plh0];
plh0=reshape(plh0,[3 numPoints])';

daterange=[+Inf ;-Inf];
for k=1:numPoints
  daterange(1)=min(daterange(1),min(ts(k).epoch));
  daterange(2)=max(daterange(2),max(ts(k).epoch));
end
numEpochs=ceil(daterange(2))-floor(daterange(1));
epochDate=daterange(1):1:daterange(2)+.5;
epochDyears=date2dyear(epochDate);

obsData=nan(numPoints,numEpochs,3);
for k=1:numPoints
   neu=ts(k).neu;
   iepoch=round(ts(k).epoch-daterange(1))+1;
   obsData(k,iepoch,1)=neu(:,1);
   obsData(k,iepoch,2)=neu(:,2);
   obsData(k,iepoch,3)=neu(:,3);
end

%% Build the output space time matrix 

[projectId,datasetId,~]=fileparts(outputfilename);
if numel(opt.projectId) > 1
   projectId=opt.projectId;
end

st = stm(projectId,datasetId,'gnss');    

techniqueAttrib=[];
techniqueAttrib.system='GNSS';
techniqueAttrib.mode='CORS';
st.techniqueAttrib=techniqueAttrib;

datasetAttrib=st.datasetAttrib;
datasetAttrib.softwareOptions=opt;
st.datasetAttrib=datasetAttrib;

st.numPoints=numPoints;
st.numEpochs=numEpochs;
st.pntName=sites(:);
st.pntCrd=[ plh0(:,1:2)*180/pi plh0(:,3) ];
st.epochDyear=epochDyears;

pntAttrib=[];
pntAttrib.events=pntEvents;
pntAttrib.dyearRange=pntDyearRange;
st.pntAttrib=pntAttrib;

epochAttrib=[];
st.epochAttrib=epochAttrib;


st.obsTypes = { 'North' 'East' 'Up' };
st.obsData = obsData.*1000;   % Convert m to mm
st.sensitivityMatrix=nan(numPoints,3,3);
st.sensitivityMatrix(:,:,1) = [ ones(numPoints,1) zeros(numPoints,1) zeros(numPoints,1) ];
st.sensitivityMatrix(:,:,2) = [ zeros(numPoints,1) ones(numPoints,1) zeros(numPoints,1) ];
st.sensitivityMatrix(:,:,3) = [ zeros(numPoints,1) zeros(numPoints,1) ones(numPoints,1) ];     
st.stochModel{1}={ ...  
            'iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sigma=0.1481,fs=24,,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-2,alphaL=-0.8281,gamma=0.992,sigma=0.06,fs=24,ndays=1,leadup=180)', ...
            'spatialcov(rho=0.0827)'};                       % GNSS North displacement
st.stochModel{2}={ ...  
            'iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sigma=0.1481,fs=24,,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-2,alphaL=-0.8281,gamma=0.992,sigma=0.06,fs=24,ndays=1,leadup=180)', ...
            'spatialcov(rho=0.1291)'};                       % GNSS East displacement
st.stochModel{3}={ ...  
            'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=1,leadup=180)' ...
            'spatialcov(rho=0.0887)'};                       % GNSS Up displacement
st.stochData=[];

% Structure array with data on each input dataset

inputDatasets=struct([]);
for k=1:numel(ts)
   inputDatasets(k).datasetId=ts(k).station;
   inputDatasets(k).techniqueId='gnss';
   datasetAttrib=[];
   datasetAttrib.softwareName= softwareName;
   datasetAttrib.softwareOptions= [];
   datasetAttrib.fileFormat=opt.inputFormat;
   datasetAttrib.fileName=ts(k).filename;
   datasetAttrib.projectId='';
   inputDatasets(k).datasetAttrib=datasetAttrib;
   inputDatasets(k).numPoints=1;
   inputDatasets(k).numEpochs=numel(ts(k).year);
   inputDatasets(k).inputDatasets=[];
end

st.inputDatasets=inputDatasets;

% Global Attributes

st.globalAttrib=opt.globalAttrib;

%% Add meteo data to space time matrix dataset

if ~isempty(opt.meteoProvider)

   meteoStations = opt.meteoStations;
   meteoBaseURL = opt.meteoBaseURL;
   meteoFileTemplate = opt.meteoFileTemplate;

   % Initialize the epoch and point attributes with meteo data
   % - The temperature and pressure for each stations are stored as epoch attributes
   % - The distance to each meteo station, which is different for each GPS
   %   station, is stored as point attribute

   numMeteoStations=numel(meteoStations);

   epochDate=dyear2date(st.epochDyear);
   
   techniqueAttrib=st.techniqueAttrib;
   epochAttrib=st.epochAttrib;
   pntAttrib=st.pntAttrib;

   epochAttrib.temperature=nan(numMeteoStations,numEpochs);
   epochAttrib.pressure=nan(numMeteoStations,numEpochs);

   pntAttrib.meteoStationDistance=zeros(numPoints,numMeteoStations);

   % Retrieve the meteo data (provide dependent)
    
   switch opt.meteoProvider
      case 'KNMI'
        % For the decomposition of the GPS signal we require temperature and air
        % pressure data. This data is for instance available at the KNMI website
        % <http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/uurgeg_280_2011-2020.zip>
        %
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
        %
        % Structure with information on meteo stations (example)
        %
        %         % Name KNMI station     Number   Coordinates    
        %         meteoStations = cell2struct({ ...
        %              'Eelde'                 280     [ 53+08/60  6+35/60 0 ] ; ... 
        %              'Nieuw Beerta'          286     [ 53+12/60  7+09/60 0 ] ; ...
        %              'Leeuwarden'            270     [ 53+13/60  5+46/60 0 ] ; ...
        %            %  'Hoorn (Terschelling)'  251     [ 53+23/60  5+21/60 0 ] ; ...
        %            %  'Stavoren'              267     [ 52+53/60  5+23/60 0 ] ; ...
        %            %  'Marknesse'             273     [ 52+42/60  5+53/60 0 ] ; ... 
        %            %  'Lauwersoog'            277     [ 53+25/60  6+12/60 0 ] ; ...
        %            %  'Hoogeveen'             279     [ 52+44/60  6+31/60 0 ] ; ...
        %            },{'name', 'number', 'coordinates' },2) ;
        %         meteoBaseURL = 'http://cdn.knmi.nl/knmi/map/page/klimatologie/gegevens/uurgegevens/';
        %         meteoFileTemplate = { 'uurgeg_###_2001-2010.txt' ; 'uurgeg_###_2011-2020.txt' ; 'uurgeg_###_2021-2030.txt' };
        
        % Download meteo data into directory opt.meteodir

        if ~isfolder(opt.meteodir)
           mkdir(opt.meteodir); 
        end
        for k=1:numMeteoStations
           for l=1:numel(meteoFileTemplate)
              meteofile=strrep(meteoFileTemplate{l},'###',num2str(meteoStations(k).number));
              if ~strcmp(meteoBaseURL,'localhost')
                 fprintf('Retrieve and unzip %s \n',meteofile)
                 urlmeteo=[ meteoBaseURL strrep(meteofile,'.txt','.zip') ];   
                 unzip(urlmeteo,opt.meteodir);
              else
                 fprintf('Meteo file %s expected in meteo data folder, make sure it is up to date!\n',meteofile)
              end
           end
        end

        % Read meteo data from file using the function |getmeteo| and save the
        % meteo data as epoch attributes. |getmeteo| optionally  plots the
        % temperature and air pressure.

        for k=1:numMeteoStations
           % read the meteo files
           meteostationname=meteoStations(k).name;
           meteofiles=cellfun(@(x) [ strrep(x,'###',num2str(meteoStations(k).number)) ],meteoFileTemplate,'UniformOutput',false);
           meteofiles=cellfun(@(x) [ opt.meteodir '/' x ],meteofiles,'UniformOutput',false);
           [meteoDate,meteoTemp,meteoPres]=getmeteo(meteostationname,meteofiles,'doplot',opt.plotMeteo,'saveplot',opt.savePlots);
           % find common epochs 
           [~,ia,ib]=intersect(floor(epochDate),round(meteoDate));
           % save in the epoch attributes
           epochAttrib.temperature(k,ia)=meteoTemp(ib);   
           epochAttrib.pressure(k,ia)=meteoPres(ib);   
        end

    otherwise
        fprintf('Unsupported meteo data provider %s, no meteo data added\n',opt.meteoProvider)
        numMeteoStations=0;
  end

  % Add meteo station names to technique attributes
  
  techniqueAttrib.meteoStations=meteoStations;

  % Distance to each meteo station
  
  for k=1:numMeteoStations
     neu = plh2neusp(st.pntCrd,meteoStations(k).coordinates);
     pntAttrib.meteoStationDistance(:,k)=sqrt(neu(:,1).^2+neu(:,2).^2);
  end
  
  % Save attributes to space time matrix
  
  st.techniqueAttrib=techniqueAttrib;
  st.epochAttrib=epochAttrib;
  st.pntAttrib=pntAttrib;

end

%% Save the integrated dataset

stmwrite(st,outputfilename);

% Finish the function

fprintf('%s finished at %s  (elapsed time %.2f s)\n',progname,datestr(now),toc)
diary OFF

catch ME

   getReport(ME)
   fprintf('%s ABORTED with an error at %s (elapsed time %.2f s)\n',progname,datestr(now),toc)
   diary OFF
   rethrow(ME)
    
end

end

%% Local functions (have their own workspace)

function [nodata,interv,intervdate,nodaydata]=namnodata(epoch,neu,crit,maxhour)
%NAMNODATA   Check timeseries for periods without new data.
%   NODATA=NAMNODATA(EPOCH,NEU) detects days in the NEU timeseries
%   that the North and East data does not appear to change. This is
%   a byproduct of the Filter approach used by SSRPOST, and this data
%   must be removed for the correct statistics.

if nargin < 1
  error('not enough arguments')
end
if ( isstruct(epoch) || ischar(epoch) )
  if isstruct(epoch)
     tseries=epoch;
  else
     load([ epoch '.mat']);
  end
  station=tseries.station;
  epoch=tseries.epoch;
  neu=tseries.neu;
%elseif nargin < 3
else
   station='';
%else
%  error('not enough arguments')
end
if nargin < 3
  crit=.0002;
end
if nargin < 4
  maxhour=18;
end

% compute change between data points (using lat and lon only)

dneu=diff(neu);
%d=sqrt(dneu(:,1).^2+dneu(:,2).^2+dneu(:,3).^2);
d2=sqrt(dneu(:,1).^2+dneu(:,2).^2);
id2=d2 < crit;

% compute start and stop of intervals with missing data 
%    start 0 -> 1
%    end 1 -> 0

idd2=[0 ; diff([0; id2]) ];

interv1=find(idd2 == 1);
interv2=find(idd2 == -1)-1;

l1=size(interv1,1);
l2=size(interv2,1);
if l2 < l1
  interv2=[interv2;size(idd2,1)];
elseif l1 > l2
  error(['l1 > l2 : ' num2str(l1) ' > ' num2str(l2) ]);
end

% remove "false alarms"...

numinterv=interv2-interv1+1;
interv1(numinterv < maxhour)=[];
interv2(numinterv < maxhour)=[];
numinterv(numinterv < maxhour)=[];

% merge intervals (i.e. eliminate short data periods

kdel=0;
for k=2:size(numinterv)
   kk=k-kdel;
   if interv1(kk)-interv2(kk-1) <= maxhour
     interv2(kk-1)=interv2(kk);
     interv1(kk)=[];
     interv2(kk)=[];
     kdel=kdel+1;
   end     
end
numinterv=interv2-interv1+1;

interv=[ interv1 interv2 numinterv ];
intervdate=[ epoch(interv1) epoch(interv2) epoch(interv2)-epoch(interv1)+1/24 ];
nodata=false(size(epoch));
for k=1:size(interv,1)
  nodata(interv(k,1):interv(k,2))=true;
end

% print message
for k=1:size(interv,1)
  fprintf('removed %s - %s  (%d epochs)\n',datestr(epoch(interv(k,1)),0),datestr(epoch(interv(k,2)),0),interv(k,2)-interv(k,1)+1)
end

% compute continuous day of year number and days with missing data

[year0,~,~]=datevec(epoch(1));
cdoy=floor(epoch-datenum(year0,1,1))+1;
[mdoy,sdoy]=grpstats([0 ;id2],cdoy,{'mean','gname'});
sdoy=str2num(char(sdoy));
nodaydata=sdoy(mdoy > .95);

% if not output arguments, plot the data

if nargout < 1

  interv
  
  figure;
  plot(epoch(2:end),[dneu(:,1)-.001 dneu(:,2) dneu(:,3)+.001]);
  title([ 'First differences ' station])
  datetick('x','mmm')
  ylabel('\Delta [m]')
  hold on
  plot(epoch(find(nodata)),nodata(find(nodata))*.0005,'kx')
  plot(epoch(find(nodata)),nodata(find(nodata))*-.0005,'kx')
  legend('N','E','U','nodata')

end

end




