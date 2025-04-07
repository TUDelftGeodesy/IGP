function varargout=stmdecomposegnss(inputfilename,outputfilename,varargin)
%stmdecomposegnss   Decompose GNSS time series.
%   STMDECOMPOSEGNSS(INPUTFILENAME,OUTPUTFILENAME,OPTIONS) does a time series
%   decomposition of the GNSS space time matrix dataset INPUTFILENAME and
%   creates an output space time matrix OUTPUTFILENAME with the decomposion
%   using the specifications given in OPTIONS. OPTIONS is a cell 
%   array or structure with the processing options, if empty, or when 
%   fields are missing, default values are used.
%
%   STMDECOMPOSEGNSS(...,UPDATEFLAG) affects the processing in case OUTPUTFILENAME
%   exists. The following values of UPDATEFLAG are accepted
%
%     'create'     abort processing if the output dataset exists
%     'update'     only do processing if one or more input dataset are newer
%                  than the output dataset
%     'overwrite'  overwrite the output dataset if it already exists
%
%   STAT=STMDECOMPOSEGNSS(...) returns a status code STAT. A status code of 0
%   indicates success.
%
%   Options for the decomposition can be given in the structure, or cell 
%   array with option/value pairs, OPTIONS. Supported options are
%
%      verbose           Verbosity level, higher is more output, 0 is almost nothing (default 0)
%      doplots           Plot level, 0 is no plots, higher is more detailed plotting (default 0)
%
%      inputDir          Directory with the inputfiles (Default '')
%
%      ROI=[]            Region of interest, as [latmin lonmin ; latmax lonmax] bounding 
%                        box, or lat/lon polygon, or kml/shape file (Default all)
%      POI=[-Inf +Inf]   Period of interest [ dYearStart dYearEnd ] or 
%                        [ dYearStart dYearEnd ; dYearStart dYearEnd ;...  ]  (Default all)
%
%      stationEvents     Name of the file with the station receiver name, antenna name, steps and 
%                        other events, only needed in case the events in the input stm needs to
%                        be overwritten (default none)
%      maxresid          Method for outlier detection, or maximum value for the residuals [N E U]  
%                        (default 'rmafilt')
%
%      usemeteo          Use meteo data (if available) in decomposition (default true), the meteo data 
%                        itself should be in the input stm
%      harmonics         Harmonic components in decomposition (default [1 1/2])
%      minyearharmonic   Minimum length of series (in year) to estimate harmonic components (default 1)
%
%      excludeFromCM     Cell array with station names to exclude from common mode and periodogram 
%                        computations (default empty)
%      removeCM          If true, do a second iteration to remove the common mode effects (default false)
%
%      includeHarmonics  Cell array with station names for which to include harmonic components in the 
%                        output space time matrix (default empty)
%      rmafilt           Number of days to be used for robust moving average filter in output space
%                        time matrix (default 0), if less or equal to one, the moving average filter 
%                        is not applied
%      datasetId         DatasetId for the output data space time matrix, if empty, "_reduced" is added 
%                        to the datasetId of the input space time matrix
%
%      saveintermediate  If true, save intermediate results in mat file (default false)
%      saveplotdir       If not empty, save plots as png in this folder (default empty) 
%
%      globalAttrib      Struct with global attribute updates (default from input file)
%      projectId         Project Id (default directory of the outputfile, or empty if no directory is specified)
%
%   Station events are already included in the input space time matrix. However, if you want to override 
%   these values you can specify an optional station events (stationEvents) file.
%
%   The default outlier detection method is a robust running moving average filter (using medians). 
%   Outliers are removed automatically, and are reported to the log file. Possible steps are only 
%   identified, it is up to the user to add these to the station events file.
%
%   If temperature influence and atmospheric loading is to be estimated, then the temperature and 
%   pressure data must be included as epoch atrributes in the space time matrix.
% 
%   Examples:
%      stmdecomposegnss('groningen_GNSS.mat','groningen_GNSS_decomposed.mat')      
%      stmdecomposegnss('groningen_GNSS.mat','groningen_GNSS_decomposed.mat',options,'update')      
%
%   See also GNSS2STM, STMREDUCEGNSS, STM, STMCHECKARGUMENTS, STMREAD and STMDIFF.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:   5 November 2020 by Hans van der Marel
% Modified: 10 March 2021 by Hans van der Marel 
%              - first alpha release
%              - working version for two examples and options in igpReduce06gps
%              - systematic testing of all possible options pending
%           29 August 2021 by Hans van der Marel
%              - Split the decomposition and reduction parts in two
%                functions
%           21 September 2021 by Hans van der Marel
%              - Fixed bug in datasetAttrib written to stm
%              - Added options to datasetAttrib
%              - Changed behaviour of rma (only trend+residual)
%              - Added includeHarmonic and rmafilt to pntAttrib
%              - residuals stack plotting under control of doplot option
%                (modified tseries toolbox function)
%              - remove pntId attribute warning
%           24 October 2021 by Hans van der Marel
%              - use getpntmask and getepochmask (removed local function bbox2poly) 
%              - default station exclusion is now empty

% Known issues: 
%              - common mode removal only works in case meteo data is
%                present, removeCM=true will not work (due to limitations 
%                of tseries toolbox)
%              - this function is a wrapper around the tseries toolbox,
%                using a dedicated structure, converting between the stm
%                and the tseries structure, making it more complicated,
%                rewrite needed to be able to use the stm sructure
%                directly, which would result in much cleaner and
%                elegant code  

%% Set constants

[a,f] = inqell('WGS-84');
e2 = 2*f - f^2;

%% Check the input arguments and options

progname='stmdecomposegnss';

%if nargin < 2
%   error('This function expects at least two input arguments.')
%end

% Default options

opt.verbose=0;                  % Default verbosity level, higher is more output, 0 is almost nothing
opt.doplots=0;                  % Default plot level, 0 is no plots, higher is more detailed plotting

opt.inputDir='';                % Directory with the inputfiles
opt.stationEvents='';           % Name of the file with the station receiver name, antenna name, steps and other events

opt.ROI=[];                     % Region of interest, as [latmin lonmin ; latmax lonmax] bounding box, or lat/lon polygon, kml/shape file (Default none)
opt.POI=[-Inf +Inf];            % Period of interest [ dYearStart dYearEnd ] or [ dYearStart dYearEnd ; dYearStart dYearEnd ; ... ] (Default none)

opt.usemeteo=true;              % Use meteo data (if available) in decomposition
opt.harmonics=[1 1/2];          % Harmonic components in decomposition (default [1 1/2])
opt.minyearharmonic=1;          % Minimum length of series (in year) to estimate harmonic components

opt.maxresid='rmafilt';         % Method for outlier detection, or maximum value for the residuals [N E U]  (default 'rmafilt')

opt.excludeFromCM={};           % Stations to exclude from common mode and periodogram computations
opt.removeCM=false;             % If true, do a second iteration to remove the common mode effects (default false)

opt.includeHarmonics={};        % Stations to include harmonic components in final estimate
opt.rmafilt=0;                  % Number of days for final robust moving average filter (default 0), if less or equal to one, no filtering

opt.datasetId='';               % DatasetId for the output data set, if empty, "_reduced" is added to the input datasetId.
opt.projectId='';               % Project Id (usualy empty when not part of integration project)

opt.saveintermediate=false;     % If true, save intermediate results in mat file (default false)
opt.saveplotdir='';             % If not empty, save plots as png in this folder (default empty) 

opt.globalAttrib=struct([]);    % Struct with global attribute updates (default from input file)


%% Duplicate output to file, and catch errors, and start timing

try

[~,outputfileroot]=fileparts(outputfilename);
diary([ outputfileroot '_' datestr(now,30) '.log' ])

fprintf('%s started at %s\n',progname,datestr(now));
tic;

% Check the options and if necessary overwrite the default values
%
% - on output opt contains the merged default an changed option values
% - on output outputfilename is empty in case no processing has to be done
%   (based on updateflag)
% - in case you have defined multiple input arguments, merge them into
%   a cell array for validation (make it a cell anyway for the checking
%   to work, if you pass inputfilenames as string, the function will
%   think it is a file with input names,

[inputfilename,outputfilename,opt]= ...
    stmcheckarguments({inputfilename},outputfilename,opt,varargin{:});
if isempty(outputfilename)
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    error('A name for the outputfile is expected.');
end
if numel(inputfilename) ~= 1
    fprintf('%s aborted abnormally at %s\n',progname,datestr(now));
    error('Only one inputfile is expected.');
end
inputfilename=char(inputfilename);

%% Get projectId 

[projectId,datasetId,~] = fileparts(outputfilename);
if numel(opt.projectId) > 1
   projectId=opt.projectId;
end

%% Read the GNSS STM

fprintf('\nReading GNSS STM ... \n');

st = stmread(fullfile(opt.inputDir,char(inputfilename)));

%% Set some frequently used variables and optionally select subset of points and epochs

numPoints=st.numPoints;
numEpochs=st.numEpochs;

pntCrd=st.pntCrd;
pntAttrib=st.pntAttrib;
pntEvents=pntAttrib.events;

epochDyears=st.epochDyear;
epochDate=dyear2date(epochDyears);
epochAttrib=st.epochAttrib;

% Select region(s) of interest (ROI) and period(s) of interest (POI)

pntMask = getpntmask(pntCrd,opt.ROI);
epochMask = getepochmask(epochDyears,opt.POI);


%% Prepare meteo data
%
% First, find out if meteo data has been stored in the space time matrix,
% if so the following fields are expected in the point and epoch attributes
%
%     epochAttrib.temperature
%     epochAttrib.pressure
%     pntAttrib.meteoStationDistance

havetemp=false;
havepressure=false;
numMeteoStations=0;
if isfield(epochAttrib,'temperature')
   havetemp=true;
   numMeteoStations=size(epochAttrib.temperature,1);
end
if isfield(epochAttrib,'pressure')
   havepressure=true;
   numMeteoStations=size(epochAttrib.pressure,1);
end

if isfield(pntAttrib,'meteoStationDistance') 
   if size(pntAttrib.meteoStationDistance,2) ~= numMeteoStations
      error('Number of meteo stations in point and epoch attributes do not match.')
   end
elseif havetemp || havepressure
   pntAttrib.meteoStationDistance=ones(numPoints,numMeteoStations);
end

if numMeteoStations <= 0 && opt.usemeteo
   fprintf('Option to use meteo data was set, but no meteo data available from input space time matrix dataset.')
   opt.usemeteo=false;
end

if opt.usemeteo && ( ~havetemp || ~havepressure )
   error('Sorry, you need both temperature and pressure data in this version of the software')
end


%% Convert space time matrix dataset into a tseries structure array ts

kk=0;
for k=1:numPoints

    % Skip points not in ROI

    if ~pntMask(k), continue; end
    
    % Store time series data
    
    kk=kk+1;

    ts(kk).station=st.pntName{k};                % station name
    plh0=[ pntCrd(k,1:2)*pi/180 pntCrd(k,3) ];
    ts(kk).plh0=plh0;                            % ref. latitude, longitude and height ([rad] and [m])
    ts(kk).vneu0=zeros(1,3);                     % a-priori velocity [m/y]

    iepoch = epochMask & all(~isnan(st.obsData(k,:,:)),3);  % epochs with data 
       
    ts(kk).t0=round(mean(epochDyears(iepoch)));  % reference time in years (must be rounded to whole years!!)
    ts(kk).epoch=epochDate(iepoch)';             % epoch as Matlab datenumber
    ts(kk).year=epochDyears(iepoch)';            % epoch in decimal years 

    events=struct([]);                           % structure array with events 
    ctmp=strsplit(pntEvents{k},';');   
    for l=1:numel(ctmp)
       ctmp2=strsplit(ctmp{l},',');
       events(l).type = ctmp2{1};
       events(l).year=date2dyear(datenum(ctmp2{2},'yyyymmddTHHMM'));
       events(l).name = ctmp2{3};
    end
    ts(kk).events=events;

    N = a ./ sqrt(1 - e2 .* sin(plh0(1)).^2);    % Conversion factors NEU <-> PLH
    M = N * (1 -e2) / ( 1 - e2 .* sin(plh0(1)).^2 );
    flat= M;
    flon= N*cos(plh0(1));

    neu= [ st.obsData(k,iepoch,1)' st.obsData(k,iepoch,2)' st.obsData(k,iepoch,3)' ]./1000;
    plh = [ plh0(1)+neu(:,1)/flat  plh0(2)+neu(:,2)/flon  plh0(3)+neu(:,3) ];
 
    ts(kk).neu=neu;                              % neu time series
    ts(kk).sneu=0.001*ones(size(neu));           % a-priori standard deviation

    ts(kk).geo2m=[flat flon 1];                  % conversion factors for geographic coordinates [deg] to meters
    ts(kk).plh=plh;

    ts(kk).filename=inputfilename;

    % Store meteo data, interpolate from data in space time matrix 

    if numMeteoStations > 0
       % Compute meteo station weights for interpolation (inverse distance)
       mwght=1./pntAttrib.meteoStationDistance(k,:);
       % Interpolate the temperature data
       if havetemp
          msel=all(isfinite(epochAttrib.temperature),2);
          ts(kk).temperature=epochAttrib.temperature(msel,iepoch)'*mwght(msel)'/sum(mwght(msel));
       end
       % Interpolate the pressure data
       if havepressure
          msel=all(isfinite(epochAttrib.pressure),2);
          ts(kk).pressure=epochAttrib.pressure(msel,iepoch)'*mwght(msel)'/sum(mwght(msel));
       end
    end
end

%% Overwrite events with events from event file (optional)
%
% - events are already stored in the input space time matrix (and ts), but may not be up to date,
%   especially if new jumps are detected
% - current events can be overwritten by events from a (updated) events file

if ~isempty(opt.stationEvents)
   eventfile=fullfile(opt.inputDir,opt.stationEvents);
   if exist(eventfile,'file')
      for k=1:numel(ts)
         ts(k).events=getEvents(eventfile,ts(k).station,ts(k).year(1),Inf);
      end
   else
      error('Incorrect stationEvent file, please specify the correct name.')     
   end
end

%% GPS decomposition (fitting) - First iteration
% The GPS timeseries is decomposed into several components: _trend_ , 
% _Atmospheric loading_, _Temperature Influence_, _Periodic components_ and 
% _Residuals_ (unmodelled effects). These components are estimated using
% the function |tseriesanalysis|, which is called from the function
% |tseriesfit|. The results are written to a timeseries structure array |tsfit|  
% for each call of |tseriesfit|. 
%
% The trendmodel is selected automatically in |tseriesfit|. For time series
% shorter than two years a linear fit is chosen, for time series over two
% year a spline fit is used. Meteo data from Eelde is used to estimate
% the temperature influence and loading effects. Also any data points
% with residuals larger than the limits set in |maxresid| are removed.

id=st.datasetId;
solid='it1';

clear tsfit;
fprintf('\nTime series decomposition %s (%s)\n',id,solid)
for k=1:numel(ts)
   fprintf('\n>>> Decompose %s timeseries for station %s:\n\n',id,ts(k).station)
   if ( ts(k).year(end)-ts(k).year(1) ) > opt.minyearharmonic && ~ isempty(opt.harmonics)
      if opt.usemeteo
         metstruct.YYYYMMDD=datestr(ts(k).epoch,'yyyymmdd');
         metstruct.Pday=ts(k).pressure;
         metstruct.Tday=ts(k).temperature;
         tsfit(k)=tseriesfit(ts(k),'harmonic',opt.harmonics,'breaks',{'ANT','DPL'},'meteo',metstruct,'maxresid',opt.maxresid);
      else
         tsfit(k)=tseriesfit(ts(k),'harmonic',opt.harmonics,'breaks',{'ANT','DPL'},'maxresid',opt.maxresid);
      end
   else
      tsfit(k)=tseriesfit(ts(k),'harmonic',[],'maxresid',opt.maxresid);
   end
end
fprintf('\n<<< Timeseries decomposition done.\n\n')

if opt.saveintermediate
   save(['tsfit' id '_' solid '.mat'],'tsfit');
end

%% Print parameter summary
% The function |tseriessummary| prints the values of the estimated 
% parameters and their standard deviation, as well as results form the 
% overall model test

if opt.doplots >= 1
   hbar=tseriessummary(tsfit,solid,id);
   if ~isempty(opt.saveplotdir)   
      plotdir=opt.saveplotdir;
      mkdir(fullfile(plotdir,[id '_' solid],'stack'));
      savepng(hbar,fullfile(plotdir,[id '_' solid],'stack'));
   end
else
   tseriessummary(tsfit,solid)
end

prtSteps(tsfit,solid)

%% Plot maps
% The function |tseriesplotmap| plots the location of the stations together
% with a representation of the estimated parameters in map format.

if opt.doplots >= 1
    
   if ~isempty(opt.harmonics)
      hmap(1)=tseriesplotmap(tsfit,'Annual',solid,id);
      %tseriesplotmap(tsfit,'AnnualP',solid,projectname);
      %tseriesplotmap(tsfit,'AnnualA',solid,projectname);
   end
   if opt.usemeteo
      hmap(2)=tseriesplotmap(tsfit,'Tempi',solid,id);
   end
   hmap(3)=tseriesplotmap(tsfit,'Velocity',solid,id);
   hmap(4)=tseriesplotmap(tsfit,'Cov',solid,id);

   if ~isempty(opt.saveplotdir)  
      plotdir=opt.saveplotdir;
      mkdir(fullfile(plotdir,[id '_' solid],'series'));
      savepng(hmap,fullfile(plotdir,[id '_' solid],'stack'));
   end

end

%% Select stations for periodogram, residual stack and common mode computation
% We will only select good stations, this means excuding not well behaving 
% stations e.g. AME2, AWG1, NORG and VEEN, given by opt.excludeFromCM

isel=find(~ismember({tsfit.station},opt.excludeFromCM));

%% GPS Periodogram - First iteration
% The Matlab function |tseriesperiodogram.m| computes the Lomb-Scargle 
% periodogram for several components. The periodogram can be computed
% both for individual stations as well as all stations together (stacked 
% periodogram). In this script we compute two periodograms: one of the 
% detrended signal, and one for the residuals. This may take some time.

if opt.doplots >= 2

   %tseriesperiodogram(stations,'tempi',solid);
   %tseriesperiodogram(tsfit,{'harmonic','tempi','atmld','residual'},solid);
   tseriesperiodogram(tsfit(isel),{'harmonic','tempi','atmld','residual'},solid);
   set(gcf,'Name',[ id '_Periodogram_Stacked_' solid ],'NumberTitle','off')
   hper(1)=gcf;

   %tseriesperiodogram(tsfit,'residual',solid);
   tseriesperiodogram(tsfit(isel),'residual',solid);
   set(gcf,'Name',[ id '_Periodogram_Residuals_' solid ],'NumberTitle','off')
   hper(2)=gcf;
   
   if ~isempty(opt.saveplotdir)  
      plotdir=opt.saveplotdir;
      mkdir(fullfile(plotdir,[id '_' solid],'stack'));
      savepng(hper,fullfile(plotdir,[id '_' solid],'stack'));
   end

end

%% Common mode - Residual stack (1) and common mode of parameters (2)
%
% The function |tseriesresidualstack| computes and plots a residual stack
% from the timeseries. The residual stack is saved to a mat file.
%
% The function |tseriescmfit| computes the common mode of the estimated
% harmonic, temperature influence and loading parameters. The common
% mode parameters are plotted using the function |tseriescmeval|.
%
% The function |tseriescmstack| is used to produce a timeseries structure
% with the common mode for plotting .

% Residual stack (the function does some plotting...)

rstack=tseriesresidualstack(tsfit(isel),solid,id,opt.doplots,opt.saveplotdir);

% Common mode of parameters

cm=tseriescmfit(tsfit(isel),solid);

% Convert residual stack and common mode into a timeseries object

if opt.usemeteo
   % We use meteo data from the first meteo station for teh cm timeseries   
   metstruct.YYYYMMDD=datestr(epochDate,'yyyymmdd');
   metstruct.Pday=epochAttrib.pressure(1,:)';
   metstruct.Tday=epochAttrib.temperature(1,:)';   
   tscm=tseriescmstack(rstack,cm,metstruct);
else
   tscm=tseriescmstack(rstack,cm);
end

% Plot common mode (w/o residual stack

if opt.doplots >= 1
   
   tseriescmeval(cm,min(rstack.ryearday):1/365:max(rstack.ryearday));
   set(gcf,'Name',[ id '_cmode_Harmonics_' solid ],'NumberTitle','off')
   hcm(1)=gcf;

   if opt.usemeteo
      tseriescmeval(cm,min(rstack.ryearday):1/365:max(rstack.ryearday),metstruct);
      %tseriescmeval(cm,epochyears,metstruct);
      set(gcf,'Name',[ id '_cmode_Harmonics+Env_' solid ],'NumberTitle','off')
      hcm(2)=gcf;
   end

   if ~isempty(opt.saveplotdir)   
      plotdir=opt.saveplotdir;
      savepng(hstack,fullfile(plotdir,[id '_' solid],'stack'));
      savepng(hcm,fullfile(plotdir,[id '_' solid],'stack'));
   end

end

% save residual stack and common mode

if opt.saveintermediate
   save(['rstack' id '_' solid '.mat'],'rstack','cm','tscm');
end

%% Plot individual components 
%
% The results from the decomposition are also be plotted component by
% component, using a single plot for all stations. .
%
% The common mode is added to the plots

%range1=[10 10 20];
range1=10;
range2=10;
rowmargin1=[1 4];
rowmargin2=[1 1];

if opt.doplots >= 2

    hraw=tseriesplotcomp([tsfit tscm],'raw','Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);
    hnostep=tseriesplotcomp([tsfit tscm],'nostep','Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);
    hfit=tseriesplotcomp([tsfit tscm],'fit','Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);
    htrend=tseriesplotcomp([tsfit tscm],'trend','Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);

    %hjumps=tseriesplotcomp([tsfit tscm],'jumps','Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);

    hharmonic=tseriesplotcomp([tsfit tscm],'harmonic','Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);

    if opt.usemeteo
       htempi=tseriesplotcomp([tsfit tscm],'tempi','Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);
       hharmonictempi=tseriesplotcomp([tsfit tscm],{'harmonic' 'tempi'},'Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);
       hatmld=tseriesplotcomp([tsfit tscm],'atmld','Id',solid,'Project',id,'Spacing',5,'RowMargin',rowmargin2);
    end

    hresidual=tseriesplotcomp([tsfit tscm],'residual','Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);

    if ~isempty(opt.saveplotdir)   
       plotdir=opt.saveplotdir;
       mkdir(fullfile(plotdir,[id '_' solid],'series'));
       savepng([hraw hnostep hfit htrend ],fullfile(plotdir,[id '_' solid],'series'));
       savepng([hatmld htempi hharmonic hharmonictempi ],fullfile(plotdir,[id '_' solid],'series'));
       savepng(hresidual,fullfile(plotdir,[id '_' solid],'series'));
   end

end

% Plot final corrected time series and removed component
%
% Two final (most important?) plots are made
% * the final corrected time series consisting of the estimated trend and residuals
% * the harmonic, temperature influence and loading components that have
%    been removed from the original series

if opt.doplots >= 1

    htrendresidual=tseriesplotcomp([tsfit tscm],{'trend' 'residual'},'Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);
    hharmonictempiatmld=tseriesplotcomp([tsfit tscm],{'harmonic' 'tempi' 'atmld'},'Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);

    if ~isempty(opt.saveplotdir)   
       plotdir=opt.saveplotdir;
       mkdir(fullfile(plotdir,[id '_' solid],'series'));
       savepng([htrendresidual hharmonictempiatmld],fullfile(plotdir,[id '_' solid],'series'));
    end

end


%% GPS decomposition (fitting) - Second iteration 
%
% The residual stack and parameter common mode can be removed from
% the time series and the time series fit can be repeated a second time.

if opt.removeCM

    solid='it2';

    % The results from the previous iterations are assigned to |tsfit1|.

    tsfit1=tsfit;
    tscm1=tscm;

    % Do the decomposition again, with common modes removed

    clear tsfit;
    fprintf('\nTime series decomposition %s (%s)\n',id,solid)
    for k=1:numel(ts)
       fprintf('\n>>> Decompose %s timeseries for station %s:\n\n',id,ts(k).station)
       if ( ts(k).year(end)-ts(k).year(1) ) > opt.minyearharmonic && ~ isempty(opt.harmonics)
          if opt.usemeteo
             metstruct.YYYYMMDD=datestr(ts(k).epoch,'yyyymmdd');
             metstruct.Pday=ts(k).pressure;
             metstruct.Tday=ts(k).temperature;
             tsfit(k)=tseriesfit(ts(k),'harmonic',opt.harmonics,'breaks',{'ANT','DPL'},'meteo',metstruct,'rstack',rstack,'cm',cm,'maxresid',opt.maxresid);
          else
             tsfit(k)=tseriesfit(ts(k),'harmonic',opt.harmonics,'breaks',{'ANT','DPL'},'rstack',rstack,'cm',cm,'maxresid',opt.maxresid);
          end
       else
          tsfit(k)=tseriesfit(ts(k),'harmonic',[],'rstack',rstack,'cm',cm,'maxresid',opt.maxresid);
       end
    end
    fprintf('\n<<< Timeseries decomposition done.\n\n')

    if opt.saveintermediate
       save(['tsfit' id '_' solid '.mat'],'tsfit');
    end

    % print summary and bar plots

    if opt.doplots >= 1
       hbar=tseriessummary(tsfit,solid,id);
       if ~isempty(opt.saveplotdir)   
          plotdir=opt.saveplotdir;
          mkdir(fullfile(plotdir,[id '_' solid],'stack'));
          savepng(hbar,fullfile(plotdir,[id '_' solid],'stack'));
       end
    else
       tseriessummary(tsfit,solid)
    end

    % Plot maps

    if opt.doplots >= 1
       if ~isempty(opt.harmonics)
          hmap(1)=tseriesplotmap(tsfit,'Annual',solid,id);
          %tseriesplotmap(tsfit,'AnnualP',solid,projectname);
          %tseriesplotmap(tsfit,'AnnualA',solid,projectname);
       end
       if opt.usemeteo
          hmap(2)=tseriesplotmap(tsfit,'Tempi',solid,id);
       end
       hmap(3)=tseriesplotmap(tsfit,'Velocity',solid,id);
       hmap(4)=tseriesplotmap(tsfit,'Cov',solid,id);  

       if ~isempty(opt.saveplotdir)  
          plotdir=opt.saveplotdir;
          mkdir(fullfile(plotdir,[id '_' solid],'series'));
          savepng(hmap,fullfile(plotdir,[id '_' solid],'stack'));
       end
       
    end

    % GPS Periodogram

    if opt.doplots >= 2

       tseriesperiodogram(tsfit(isel),{'harmonic','tempi','atmld','residual'},solid);
       set(gcf,'Name',[ id '_Periodogram_Stacked_' solid ],'NumberTitle','off')
       hper(1)=gcf;

       tseriesperiodogram(tsfit(isel),'residual',solid);
       set(gcf,'Name',[ id '_Periodogram_Residuals_' solid ],'NumberTitle','off')
       hper(2)=gcf;
         
       if ~isempty(opt.saveplotdir)  
          plotdir=opt.saveplotdir;
          mkdir(fullfile(plotdir,[id '_' solid],'stack'));
          savepng(hper,fullfile(plotdir,[id '_' solid],'stack'));
       end

    end

    % Residual stack 

    rstack=tseriesresidualstack(tsfit(isel),solid,id,opt.doplots,opt.saveplotdir);

    % Common mode of parameters

    cm=tseriescmfit(tsfit(isel),solid);

    % Convert residual stack and common mode into a timeseries object

    if opt.usemeteo
       % We use meteo data from the first station   
       metstruct.YYYYMMDD=datestr(epochDate,'yyyymmdd');
       metstruct.Pday=epochAttrib.pressure(1,:)';
       metstruct.Tday=epochAttrib.temperature(1,:)';   
       tscm=tsstack(rstack,cm,metstruct);
    else
       tscm=tsstack(rstack,cm);
    end

    % plot common mode

    if opt.doplots >= 1

       tseriescmeval(cm,min(rstack.ryearday):1/365:max(rstack.ryearday));
       set(gcf,'Name',[ id '_cmode_Harmonics_' solid ],'NumberTitle','off')
       hcm(1)=gcf;

       if opt.usemeteo
          tseriescmeval(cm,min(rstack.ryearday):1/365:max(rstack.ryearday),metstruct);
          %tseriescmeval(cm,epochyears,metstruct);
          set(gcf,'Name',[ id '_cmode_Harmonics+Env_' solid ],'NumberTitle','off')
          hcm(2)=gcf;
       end
       
       if ~isempty(opt.saveplotdir)   
          plotdir=opt.saveplotdir;
          savepng(hstack,fullfile(plotdir,[id '_' solid],'stack'));
          savepng(hcm,fullfile(plotdir,[id '_' solid],'stack'));
       end

    end

    % save residual stack and common mode

    if opt.saveintermediate
       save(['rstack' id '_' solid '.mat'],'rstack','cm','tscm');
    end

    % Change name of the common mode correction before adding to plot

    tscm1.station='CMC';

    % Plot several of the estimated components 
    
    if opt.doplots >= 2

        hfit=tseriesplotcomp([tsfit tscm1 tscm],'fit','Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);
        htrend=tseriesplotcomp([tsfit tscm1 tscm],'trend','Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);

        hharmonic=tseriesplotcomp([tsfit tscm1 tscm],'harmonic','Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);

        if opt.usemeteo
           htempi=tseriesplotcomp([tsfit tscm1 tscm],'tempi','Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);
           hharmonictempi=tseriesplotcomp([tsfit tscm1 tscm],{'harmonic' 'tempi'},'Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);
           hatmld=tseriesplotcomp([tsfit tscm1 tscm],'atmld','Id',solid,'Project',id,'Spacing',5,'RowMargin',rowmargin2);
        end

        hresidual=tseriesplotcomp([tsfit tscm1 tscm],'residual','Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);
        
        if ~isempty(opt.saveplotdir)   
           plotdir=opt.saveplotdir;
           mkdir(fullfile(plotdir,[id '_' solid],'series'));
           savepng([hfit htrend ],fullfile(plotdir,[id '_' solid],'series'));
           savepng([hatmld htempi hharmonic hharmonictempi ],fullfile(plotdir,[id '_' solid],'series'));
           savepng(hresidual,fullfile(plotdir,[id '_' solid],'series'));
       end
        
    end
    
    % ... and the two final (most important) ones

    if opt.doplots >= 1

        htrendresidual=tseriesplotcomp([tsfit tscm1 tscm],{'trend' 'residual'},'Id',solid,'Project',id,'Spacing',range1,'RowMargin',rowmargin1);
        hharmonictempiatmld=tseriesplotcomp([tsfit tscm1 tscm],{'harmonic' 'tempi' 'atmld'},'Id',solid,'Project',id,'Spacing',range2,'RowMargin',rowmargin2);

        if ~isempty(opt.saveplotdir)   
           plotdir=opt.saveplotdir;
           mkdir(fullfile(plotdir,[id '_' solid],'series'));
           savepng([htrendresidual hharmonictempiatmld],fullfile(plotdir,[id '_' solid],'series'));
        end

    end

end

%% Make space time matrix with final displacement series
%
% The final displacement series is the sum of
% - the trend + residuals by default
% - the trend + residuals + harmonics for stations in opt.includeharmonics (e.g. NORG)
% - trend + residuals are optionaly smoothed using a moving average filter (opt.rmafiltDays > 1)
%      opt.rmafilt=21  number of days in moving average filter
%      opt.rmacrit=5   ....
%      opt.rmastepcrit=0  step criterion (0 is ignore steps)
% - the remainder (the sum of harmonic components, if not include in obsData, atmospheric loading, 
%   temperature influence, and moving average residuals) is saved as auxialiary data in the 
%   space time matrix (steps are not saved in the auxData).

iepoch=round(epochDate-epochDate(1))+1;

obsData=nan(numPoints,numEpochs,3);
auxData=nan(numPoints,numEpochs,3);
for k=1:numel(tsfit)
   station=tsfit(k).station;
   % Build array with displacement data (harmonic components not included)
   data=tsfit(k).neutrend+tsfit(k).neuresidual;
   remainder=tsfit(k).neutempi+tsfit(k).neuatmld+tsfit(k).neuharmonic;
   % Optionally filter the data 
   if opt.rmafilt > 1
      tmp=nan(size(data));
      tmp(:,1)=rmafilt(data(:,1),opt.rmafilt,'crit',opt.rmacrit,'stepcrit',opt.rmastepcrit);
      tmp(:,2)=rmafilt(data(:,2),opt.rmafilt,'crit',opt.rmacrit,'stepcrit',opt.rmastepcrit);
      tmp(:,3)=rmafilt(data(:,3),opt.rmafilt,'crit',opt.rmacrit,'stepcrit',opt.rmastepcrit);
      remainder=remainder+data-tmp;
      data=tmp;
   end
   % For selected stations, add harmonic components to the displacements
   if ismember(station,opt.includeHarmonics)
      data=data+tsfit(k).neuharmonic;
      remainder=remainder-tsfit(k).neuharmonic;
   end
   % Insert data into space time matrix   
   [~,kk]=ismember(station,st.pntName);
   kepoch=round(dyear2date(tsfit(k).year)-epochDate(1))+1;
   [~,ia,ib]=intersect(kepoch,iepoch);
   obsData(kk,ib,1)=data(ia,1)';
   obsData(kk,ib,2)=data(ia,2)';
   obsData(kk,ib,3)=data(ia,3)';
   auxData(kk,ib,1)=remainder(ia,1)';
   auxData(kk,ib,2)=remainder(ia,2)';
   auxData(kk,ib,3)=remainder(ia,3)';
end

%% Save common mode estimates and/or corrections as epoch attributes 

epochCM=nan(3,numEpochs);
kepoch=round(dyear2date(tscm.year)-epochDate(1))+1;
data=tscm.neutrend+tscm.neuresidual+tscm.neuharmonic+tscm.neutempi+tscm.neuatmld;
[~,ia,ib]=intersect(kepoch,iepoch);
epochCM(:,ib)=data(ia,:)';

if opt.removeCM
   epochCMC=nan(3,numEpochs);
   kepoch=round(dyear2date(tscm.year)-epochDate(1))+1;
   data=tscm1.neutrend+tscm1.neuresidual+tscm1.neuharmonic+tscm1.neutempi+tscm1.neuatmld;
   [~,ia,ib]=intersect(kepoch,iepoch);
   epochCMC(:,ib)=data(ia,:)';
end

%% Build the output space time matrix 

if ~isempty(opt.datasetId)
   datasetIdOut = opt.datasetId;
else
   datasetIdOut = datasetId;
end
if strcmpi(datasetIdOut,st.datasetId)
   datasetIdOut = [ st.datasetId '_decomposed' ];
end

numPoints=sum(pntMask);
numEpochs=sum(epochMask);
   
stout = stm(projectId,datasetIdOut,'gnss');    

techniqueAttrib=st.techniqueAttrib;
techniqueAttrib.status='Decomposed';
stout.techniqueAttrib=techniqueAttrib;

datasetAttrib=stout.datasetAttrib;
datasetAttrib.softwareOptions=opt;
stout.datasetAttrib=datasetAttrib;

stout.numPoints=numPoints;
stout.numEpochs=numEpochs;
stout.pntName=st.pntName(pntMask);
stout.pntCrd=st.pntCrd(pntMask,:);
stout.epochDyear=st.epochDyear(epochMask);

pntAttrib=st.pntAttrib;
if ~isempty(pntAttrib)
   pntAttribFields = fieldnames(pntAttrib);
   for k=1:numel(pntAttribFields)
      pntAttribField=pntAttribFields{k};
      tmp=pntAttrib.(pntAttribField);
      pntAttrib.(pntAttribField)=tmp(pntMask,:);
   end
end
stout.pntAttrib=pntAttrib;

epochAttrib=st.epochAttrib;
if ~isempty(epochAttrib)
   epochAttribFields = fieldnames(epochAttrib);
   for k=1:numel(epochAttribFields)
      epochAttribField=epochAttribFields{k};
      tmp=epochAttrib.(epochAttribField);
      epochAttrib.(epochAttribField)=tmp(:,epochMask);
   end
end
epochAttrib.commonModeEstimate=epochCM(:,epochMask);
if opt.removeCM
   epochAttrib.commonModeCorrection=epochCMC(:,epochMask);
end
stout.epochAttrib=epochAttrib;

stout.obsTypes = st.obsTypes;
stout.obsData = obsData(pntMask,epochMask,:).*1000;          % Convert m to mm
stout.sensitivityMatrix = st.sensitivityMatrix(pntMask,:,:);
stout.stochModel{1}={ ...  
            'iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sigma=0.1481,fs=24,,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-2,alphaL=-0.8281,gamma=0.992,sigma=0.06,fs=24,ndays=1,leadup=180)', ...
            'spatialcov(rho=0.0827)'};                       % GNSS North displacement
stout.stochModel{2}={ ...  
            'iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sigma=0.1481,fs=24,,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-2,alphaL=-0.8281,gamma=0.992,sigma=0.06,fs=24,ndays=1,leadup=180)', ...
            'spatialcov(rho=0.1291)'};                       % GNSS East displacement
stout.stochModel{3}={ ...  
            'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=1,leadup=180)' ...
            'spatialcov(rho=0.0887)'};                       % GNSS Up displacement
stout.stochData=[];

stout.auxTypes = {'remainderNorth' 'remainderEast' 'remainderUp' };
stout.auxData = auxData(pntMask,epochMask,:).*1000;          % Convert m to mm

% Structure array with data on each input dataset

inputDatasets=struct([]);
inputDatasets(1).datasetId=st.datasetId;
inputDatasets(1).techniqueId=st.techniqueId;
datasetAttrib=st.datasetAttrib;
datasetAttrib.fileName=inputfilename;
inputDatasets(1).datasetAttrib=datasetAttrib;
inputDatasets(1).numPoints=st.numPoints;
inputDatasets(1).numEpochs=st.numEpochs;
inputDatasets(1).inputDatasets=st.inputDatasets;

stout.inputDatasets=inputDatasets;

% Global attributes (copy from input dataset and optionally change).

globalAttrib=st.globalAttrib;
allfields=fieldnames(struct(opt.globalAttrib));
for k=1:numel(allfields)
  globalAttrib.(allfields{k})=opt.globalAttrib.(allfields{k});
end
stout.globalAttrib=globalAttrib;

%% Save parameter estimates as point attributes
%
% We do this after pntMask has been applied, but tsfit can still have 
% a different order of parameters

pntAttrib=stout.pntAttrib;

[~,ia,ib]=intersect({tsfit.station}',stout.pntName);

maxknots=0;
maxsteps=0;
maxharmonics=0;
maxpar=0;
for k=1:numel(ia)
   maxknots=max(maxknots,numel(tsfit(k).knots));
   maxsteps=max(maxsteps,numel(tsfit(k).tjump));
   maxharmonics=max(maxharmonics,numel(tsfit(k).harmonic));
   maxpar=max(maxpar,tsfit(k).parindex(end));
end

pntAttrib.trendmodel = { tsfit.method }';
pntAttrib.pdegree = nan(numPoints,1);
pntAttrib.pspline = nan(numPoints,1);
pntAttrib.knots = nan(numPoints,maxknots);
pntAttrib.tsteps = nan(numPoints,maxsteps);
pntAttrib.harmonics = nan(numPoints,maxharmonics);
pntAttrib.amplitudes = nan(numPoints,maxharmonics,3);

pntAttrib.parindex(ib,1:4) = reshape([ tsfit(ia).parindex ],4,[])';
pntAttrib.parestimates = nan(numPoints,maxpar,3);
pntAttrib.parcovmatrices = nan(numPoints,maxpar,maxpar,3);
pntAttrib.parlabels = cell(numPoints,maxpar);
pntAttrib.parunits = cell(numPoints,maxpar);

for k=1:numel(ia)
  if isempty(tsfit(ia(k)).pspline)     
      pntAttrib.pdegree(ib(k),1) = tsfit(ia(k)).pdegree;
  else
      pntAttrib.pspline(ib(k),1) = tsfit(ia(k)).pspline; 
  end
  pntAttrib.knots(ib(k),1:numel(tsfit(ia(k)).knots)) = tsfit(ia(k)).knots;
  pntAttrib.tsteps(ib(k),1:numel(tsfit(ia(k)).tjump)) = tsfit(ia(k)).tjump;
  pntAttrib.harmonic(ib(k),1:numel(tsfit(ia(k)).harmonic)) = tsfit(ia(k)).harmonic;

  pntAttrib.amplitudes(ib(k),1:numel(tsfit(ia(k)).harmonic),1) = tsfit(ia(k)).amplitude(1,:) * 1000;
  pntAttrib.amplitudes(ib(k),1:numel(tsfit(ia(k)).harmonic),2) = tsfit(ia(k)).amplitude(2,:) * 1000;
  pntAttrib.amplitudes(ib(k),1:numel(tsfit(ia(k)).harmonic),3) = tsfit(ia(k)).amplitude(3,:) * 1000;
  
  pntAttrib.parestimates(ib(k),1:tsfit(ia(k)).parindex(end),1) = tsfit(ia(k)).xlat * 1000;
  pntAttrib.parestimates(ib(k),1:tsfit(ia(k)).parindex(end),2) = tsfit(ia(k)).xlon * 1000;
  pntAttrib.parestimates(ib(k),1:tsfit(ia(k)).parindex(end),3) = tsfit(ia(k)).xrad * 1000;
  pntAttrib.parlabels(ib(k),1:tsfit(ia(k)).parindex(end))= tsfit(ia(k)).parlabel ; 
  pntAttrib.parunits(ib(k),1:tsfit(ia(k)).parindex(end))= tsfit(ia(k)).parunit ; 

  pntAttrib.parcovmatrices(ib(k),1:tsfit(ia(k)).parindex(end),1:tsfit(ia(k)).parindex(end),1) = tsfit(ia(k)).qxlat * 1e6;
  pntAttrib.parcovmatrices(ib(k),1:tsfit(ia(k)).parindex(end),1:tsfit(ia(k)).parindex(end),2) = tsfit(ia(k)).qxlon * 1e6;
  pntAttrib.parcovmatrices(ib(k),1:tsfit(ia(k)).parindex(end),1:tsfit(ia(k)).parindex(end),3) = tsfit(ia(k)).qxrad * 1e6;

end

pntAttrib.meanstd(ib,1:3) = reshape([ tsfit(ia).meanstd ],3,[])' * 1000;
pntAttrib.empstd(ib,1:3) = reshape([ tsfit(ia).empstd ],3,[])' * 1000;
pntAttrib.rms(ib,1:3) = reshape([ tsfit(ia).rms ],3,[])' * 1000;
pntAttrib.omt(ib,1:3) = reshape([ tsfit(ia).omt ],3,[])';

pntAttrib.vel0(ib,1:3) = reshape([ tsfit(ia).vneu0 ],3,[])' * 1000;
pntAttrib.vel(ib,1:3) = reshape([ tsfit(ia).vel ],3,[])' * 1000;
pntAttrib.svel(ib,1:3) = reshape([ tsfit(ia).svel ],3,[])' * 1000;

pntAttrib.t0(ib,1) = reshape([ tsfit(ia).t0 ],1,[])';
pntAttrib.plh0(ib,1:3) = reshape([ tsfit(ia).plh0 ],3,[])';

% Add include harmonics flag and robust moving average parameter to the attributes

pntAttrib.includeharmonic = ismember(stout.pntName,opt.includeHarmonics); 
pntAttrib.rmafilt = ones(numPoints,1) * opt.rmafilt;

stout.pntAttrib=pntAttrib;

%% Save the decomposed dataset

stmwrite(stout,outputfilename);

%% Finish the function

fprintf('%s finished at %s  (elapsed time %.2f s)\n',progname,datestr(now),toc)
diary OFF

catch ME

   getReport(ME)
   fprintf('%s ABORTED with an error at %s (elapsed time %.2f s)\n',progname,datestr(now),toc)
   diary OFF
   rethrow(ME)
    
end

%% Nested functions (share variables with calling function; only possible when used within a function)

end

%% Local functions (have their own workspace)

% [End of main]

