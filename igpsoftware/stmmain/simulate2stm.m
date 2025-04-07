function simulate2stm(simid,t,obsScenario,obsGaps,stochModel)
%simulate2stm   Simulate 2-D deformation field with geodetic observations.
%   simulate2stm(simid,t,obsScenario,obsGaps,stochModel) simulates levelling, 
%   GNSS and InSAR displacement observations. simid is a character string that 
%   is used as base for the output space time matrices. t is a vector with the 
%   epoch times in decimal years. obsScenario, obsGaps and stochModel are either 
%   character strings with standard observation scenarions, data gap options, 
%   and stochastic model specification, or a structure defining a non-standard 
%   scenario. 
%
%   The simulation scenarios are defined at the start of the function and are
%   selected by the obsScenario, obsGaps and stochModel inputs:
%
%   obsScenario:
%      'small'             10 Lev, GNSS and InSAR (Asc 2x, Dsc 1x), fully connected 
%      'mediumSimple1'     50 Lev, GNSS and InSAR (Asc 2x, Dsc 1x), fully connected, with extra 50 Lev
%      'mediumSimple2'     50 Lev and GNSS, with 100 InSAR (Asc 2x, Dsc 1x), well connected
%      'mediumUnderObs1'   10 Lev, 5 GNSS, 90 Asc, 190 Asc and 75 Dsc InSAR, sparsely connected (low redundancy)
%      'mediumUnderObs2'   Idem, but more heavily connected
%
%   obsGaps:
%      'none'      simulation with no gaps (default)
%      'standard'  For levelling remove 1 complete epoch, and for 50% of the levelling points 
%                  remove 1 to 4 random observations. For 70% of the GNSS points, remove the first
%                  1 to 6 observations. For InSAR, keep all.
%
%   stochModel:
%      'simple'    Simple uncorrelated noise of 1 mm standard deviation
%      'standard'  Standard, correlated, noise models for each technique (default)
%
%   The obsGaps and stochModel parameters have defaults ('none' and 'standard').
%   Instead of a string, these three parameters can also be structure defining 
%   non-standard scenarios. See inside the function for guidance. 
%
%   The following output datasets are created
%
%     <simid>_<datasetId{k}>.mat         space time matrix dataset with noisy observations     
%     <simid>_<datasetId{k}>_truth.mat   space time matrix dataset with truth observations
%     <simid>_truth.mat                  space time matrix dataset with truth displacements
%
%     <simid>_filenames.txt              text file with the dataset names (noisy observation data)
%     <simid>_filenames_truth.txt        text file with the dataset names (truth observation data)
%
%   with <simid> and cell array <datasetId{k}> defined in the user input data 
%   section, and k=1...numDatasets.
%
%   The function produces a number of plots showing the simulated deformation field,
%   points and observations.
%
%   Examples:
%
%     simulate2stm('simcase0',[2003:1:2019],'mediumSimple2','none','simple')
%     simulate2stm('simcase0',[2003:1:2019],'mediumSimple2')
%
%   To obtain repeatable results use "rng('default')" to (re)set the random number 
%   generator to its default value to produce the same random numbers as if you 
%   restarted MATLAB.
%
%   See also stm, stmwrite, stmstochmodel and stmintegrate.
%
%   (c) Hans van der Marel, Delft University of Technology, 2019-2020.

%% Simulate 2-D deformation field with geodetic observations
%
% *Hans van der Marel*
%
% Simulation of observations for NAM Integrated Geodetic Processing (IGP)
% project.
%
% This function simulates levelling, GNSS and InSAR Space Time Datasets with
% simulated noise, and writes these, with the truth data, to disk. The
% function does the following
% 
% * define the observation scenario and other parameters (user defined)
% * simulate a truth vertical and horizontal velocity field 
% * simulate observation points for the different techniques using the 
%   data provided in the observation scenario
% * generate the observation datasets (truth data) 
% * randomly remove some observations to simulate missing observations (optional)
%        - for levelling remove some complete epochs and random observations
%        - for GPS, remove for some points, starting epochs
%        - for InSAR, keep all (for the time being)
% * write the 'truth' observation datasets to disk
% * simulate observation noise
% * write the observation datasets (with noise included) to disk
% * write the 'truth' displacement field to disk
%
% The following output datasets are created
%
%     <simid>_<datasetId{k}>.mat         space time matrix dataset with noisy observations     
%     <simid>_<datasetId{k}>_truth.mat   space time matrix dataset with truth observations
%     <simid>_truth.mat                  space time matrix dataset with truth displacements
%
%     <simid>_filenames.txt              text file with the dataset names (noisy observation data)
%     <simid>_filenames_truth.txt        text file with the dataset names (truth observation data)
%
%   with <simid> and cell array <datasetId{k}> defined in the user input data 
%   section, and k=1...numDatasets.
%
% The simulation scenarios are defined at the start of the function and
% can be selected based on the function arguments.
%
% Created:  15 July 2019 by Hans van der Marel (ispSimDefo2)
% Modified: 13 August 2020 by Hans van der Marel (ispSimDefo2b)
%              - output of individual datasets
%              - define multiple cases
%              - new dataset structure
%           30 September 2020 by Hans van der Marel (simulate2stm)
%              - renamed to simulate2stm
%              - convert coordinates to latitude/longitude
%              - add stochastic models and simulate errors
%              - simulate missing data
%              - create a-priory displacement dataset
%              - save truth output dataset(s)
%            3 October 2020 by Hans van der Marel (simulate2stm)
%              - moved all user defined input to the top of the script
%              - improved documentation and several minor changes
%            7 October 2020 by Hans van der Marel
%              - reworked into a function
%              - reorganized the definition section
%           29 October 2020 by Hans van der Marel
%              - simulate different start and stop times for the datasets

% Add stmutil toolbox (with stmstochmodel) to Matlab path

if isfolder('./stmutil')
   addpath('./stmutil');
elseif isfolder('../0_toolboxes/stmutil')
   addpath('../0_toolboxes/stmutil');
end

%% Check the input arguments

if nargin < 3 || nargin > 5
   error('This function must have between 3 and 5 input arguments.')
end
if nargin < 4 || isempty(obsGaps)
   obsGaps='none';
end
if nargin < 5 || isempty(stochModel)
   stochModel='standard';
end

progname='simulate2stm';
fprintf('%s started at %s\n',progname,datestr(now));
tic;

%% Define observation scenarios, observation gaps and stochastic models  
%
% In this section the techniques and number of observation points for
% each technique are specified. This is done in the form of several
% matrices and variables defining the observation scenario. 

% Define the various observation scenarios ...
%
% Points can be observed by single or multiple techniques. Most points
% will however be observed by a single technique, but a few points,
% may be observed by multiple techniques (e.g IGRS). 
%
% For new observation scenarios, copy one of the cases, rename, and modify 
% the parameters.

if isstruct(obsScenario)
    datasetId=obsScenario.datasetId;              % Cell array with solution ID
    techniqueId=obsScenario.techniqueId;          % Cell array with the technique
    incAngle=obsScenario.incAngle;                % Incidence angle in degrees (InSAR only) 
    azAngle=obsScenario.azAngle;                  % Azimuth angle in degrees (InSAR only)
    plotSymbol=obsScenario.plotSymbol;            % Symbol for plots
    nptech=obsScenario.nptech;                    % Number of points per technique
    npdualtech=obsScenario.npdualtech;            % Number of common points between techniques (lower triangular matrix)
    npalltech=obsScenario.npalltech;              % Number of points observed by all techniques
    obsScenario='other';
end
switch obsScenario
    
    case 'small'

        % small number of points with fully connected observations

        datasetId=    { 'lev'  'gps'  'sarDsc1'  'sarAsc1' 'sarAsc2'};     % Cell array with solution ID
        techniqueId=  { 'lev' 'gnss' 'insar' 'insar' 'insar' };     % Cell array with the technique
        incAngle=     [  nan    nan    38     30     42    ];    % Incidence angle in degrees (InSAR only) 
        azAngle=      [  nan    nan    100  -100    -100   ];    % Azimuth angle in degrees (InSAR only)
        plotSymbol=   {   's'   'o'    '*'   'x'    '+'    };    % Symbol for plots

        nptech=       [  10     10     10     10     10    ];    % Number of points per technique
        npdualtech=   [   0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0      0      0      0    ];    % Number of common points between techniques (lower triangular matrix)
        npalltech=    10;                                        % Number of points observed by all techniques

        % define the measurement epochs
 
        t=[2003:1:2019];
        
    case 'mediumSimple1'

        % medium number of points, with some not observed North and East components, but otherwise fully connected observations

        datasetId=    { 'lev'  'gps'  'sarDsc1'  'sarAsc1' 'sarAsc2'};     % Cell array with solution ID
        techniqueId=  { 'lev' 'gnss' 'insar' 'insar' 'insar' };     % Cell array with the technique
        incAngle=     [  nan    nan    38     30     42    ];    % Incidence angle in degrees (InSAR only) 
        azAngle=      [  nan    nan    100  -100    -100   ];    % Azimuth angle in degrees (InSAR only)
        plotSymbol=   {   's'   'o'    '*'   'x'    '+'    };    % Symbol for plots

        nptech=       [ 100     50     50     50     50    ];    % Number of points per technique
        npdualtech=   [   0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0      0      0      0    ];    % Number of common points between techniques (lower triangular matrix)
        npalltech=     50;                                        % Number of points observed by all techniques

    case 'mediumSimple2'

        % medium number of points, with some not observed North and East components, but otherwise fully connected observations

        datasetId=    { 'lev'  'gps'  'sarDsc1'  'sarAsc1' 'sarAsc2'};     % Cell array with solution ID
        techniqueId=  { 'lev' 'gnss' 'insar' 'insar' 'insar' };     % Cell array with the technique
        incAngle=     [  nan    nan    38     30     42    ];    % Incidence angle in degrees (InSAR only) 
        azAngle=      [  nan    nan    100  -100    -100   ];    % Azimuth angle in degrees (InSAR only)
        plotSymbol=   {   's'   'o'    '*'   'x'    '+'    };    % Symbol for plots

        nptech=       [  50     50    100    100    100    ];    % Number of points per technique
        npdualtech=   [   0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0      0      0      0    ; ...
                          0      0     50      0      0    ; ...
                          0      0      0      0      0    ];    % Number of common points between techniques (lower triangular matrix)
        npalltech=     50;                                        % Number of points observed by all techniques

    case 'mediumUnderObs1'

        % simulation with very few connections between the InSAR datasets (this is the original default simulation)

        datasetId=    { 'lev'  'gps'  'sarDsc1'  'sarAsc1' 'sarAsc2'};     % Cell array with solution ID
        techniqueId=  { 'lev' 'gnss' 'insar' 'insar' 'insar' };     % Cell array with the technique
        incAngle=     [  nan    nan    38     30     42    ];    % Incidence angle in degrees (InSAR only) 
        azAngle=      [  nan    nan    100  -100    -100   ];    % Azimuth angle in degrees (InSAR only)
        plotSymbol=   {   's'   'o'    '*'   'x'    '+'    };    % Symbol for plots

        nptech=       [  10      5     90    190     75    ];    % Number of points per technique
        npdualtech=   [   0      0      0      0      0    ; ...
                          2      0      0      0      0    ; ...
                          3      2      0      0      0    ; ...
                          1      0      3      0      0    ; ...
                          1      0      4      2      0    ];    % Number of common points between techniques (lower triangular matrix)
        npalltech=    1;                                         % Number of points observed by all techniques

    case 'mediumUnderObs2'
        
        % simulation with roughly 50% connectivity between the InSAR datasets
        
        datasetId=    { 'lev'  'gps'  'sarDsc1'  'sarAsc1' 'sarAsc2'};     % Cell array with solution ID
        techniqueId=  { 'lev' 'gnss' 'insar' 'insar' 'insar' };     % Cell array with the technique
        incAngle=     [  nan    nan    38     30     42    ];    % Incidence angle in degrees (InSAR only) 
        azAngle=      [  nan    nan    100  -100    -100   ];    % Azimuth angle in degrees (InSAR only)
        plotSymbol=   {   's'   'o'    '*'   'x'    '+'    };    % Symbol for plots

        nptech=       [  10      5     90    190     75    ];    % Number of points per technique
        npdualtech=   [   0      0      0      0      0    ; ...
                          2      0      0      0      0    ; ...
                          3      2      0      0      0    ; ...
                          1      0     60      0      0    ; ...
                          1      0     60     30      0    ];    % Number of common points between techniques (lower triangular matrix)
        npalltech=    1;                                         % Number of points observed by all techniques
    
    case 'other'
    otherwise
        error('unknown simulation case, aborting...')
end

% Define the observation gaps ...

if isstruct(obsGaps)
    removeobs=true;                                   % If true, simulate missing observations, using the parameters below. If false, the parameters below are ignored
    percMissingLev=obsGaps.percMissingLev;            % Percentage of points with missing levelling observations
    maxMissingLev=obsGaps.maxMissingLev;              % Maximum number of missing observations (for levelling points with missing observations)
    numMissingEpochsLev=obsGaps.numMissingEpochsLev;  % Number of missing levelling epochs / campaigns
    percMissingGNSS=obsGaps.percMissingGNSS;          % Percentage of points with missing GNSS data at start
    maxMissingGNSS=obsGaps.maxMissingGNSS;            % Maximum number of missing GNSS datapoints at start (for GNSS points with missing data)
    obsGaps='other';
end
switch obsGaps
    case 'none'
        removeobs=false;         % If true, simulate missing observations, using the parameters below. If false, the parameters below are ignored
        
        percMissingLev=0;        % Percentage of points with missing levelling observations
        maxMissingLev=0;         % Maximum number of missing observations (for levelling points with missing observations)
        numMissingEpochsLev=0;   % Number of missing levelling epochs / campaigns

        percMissingGNSS=0;       % Percentage of points with missing GNSS data at start
        maxMissingGNSS=0;        % Maximum number of missing GNSS datapoints at start (for GNSS points with missing data)
        
        maxStartDelay=0;         % simulate a delay of at most maxStartDelay at the start of datasets
        maxEarlyEnd=0;           % simulate an early end of at most maxEarlyEnd at the end of datasets

    case 'standard'
        removeobs=true;         % If true, simulate missing observations, using the parameters below. If false, the parameters below are ignored

        percMissingLev=50;      % Percentage of points with missing levelling observations
        maxMissingLev=4;        % Maximum number of missing observations (for levelling points with missing observations)
        numMissingEpochsLev=1;  % Number of missing levelling epochs / campaigns

        percMissingGNSS=70;     % Percentage of points with missing GNSS data at start
        maxMissingGNSS=6;       % Maximum number of missing GNSS datapoints at start (for GNSS points with missing data)

        maxStartDelay=3;         % simulate a delay of at most maxStartDelay at the start of datasets
        maxEarlyEnd=3;           % simulate an early end of at most maxEarlyEnd at the end of datasets

    case 'other'
    otherwise
        error('unknown observation gap model, aborting...')
end

% Define the stochastic models...

if isstruct(stochModel)
    stochModelLev=stochModel.Lev;
    stochModelGNSS=stochModel.GNSS;
    stochModelInSAR=stochModel.InSAR;
    stochModel='other';
end
switch stochModel
    
    case 'simple'

        stochModelLev{1}={'WN(sigma=1)','spatialcov(rho=Inf)'};  
        stochModelGNSS{1}={'WN(sigma=1)','spatialcov(rho=Inf)'};
        stochModelGNSS{2}={'WN(sigma=1)','spatialcov(rho=Inf)'};
        stochModelGNSS{3}={'WN(sigma=1)','spatialcov(rho=Inf)'};
        stochModelInSAR{1}={'WN(sigma=1)','spatialcov(rho=Inf)'};

    case 'standard'

        stochModelLev{1}={'WN(sigma=1)','spatialcov(rho=.9)'};    % For levelling the model is the same for each epoch
        stochModelGNSS{1}={  ...
            'iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sigma=0.1481,fs=24,,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-2,alphaL=-0.8281,gamma=0.992,sigma=0.06,fs=24,ndays=1,leadup=180)', ...
            'spatialcov(rho=0.0827)'};                       % GNSS North displacement
        stochModelGNSS{2}={  ...
            'iGGM(alphaR=-1.91,alphaL=-0.70,gamma=0.9137,sigma=0.1481,fs=24,,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-2,alphaL=-0.8281,gamma=0.992,sigma=0.06,fs=24,ndays=1,leadup=180)', ...
            'spatialcov(rho=0.1291)'};                       % GNSS East displacement
        stochModelGNSS{3}={  ...
            'iGGM(alphaR=-1.91,alphaL=-0.74,gamma=0.9861,sigma=0.0967,fs=24,ndays=1,leadup=180)', ...
            'iGGM(alphaR=-1.91,alphaL=-0.8525,gamma=0.999096,sigma=0.04,fs=24,ndays=1,leadup=180)' ...
            'spatialcov(rho=0.0887)'};                       % GNSS Up displacement
        stochModelInSAR{1}={'tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)'};   % InSAR

    case 'other'
    otherwise
        error('unknown stochastic model')      
end

% End of definitions

%% Check consistency of nptech, npdualtech and npalltech 

ok=true;
for l=1:numel(techniqueId)
   if sum(npdualtech(:,l))+sum(npdualtech(l,:))+npalltech > nptech(l)
      warning(['Too many dual technique observation points for ' datasetId{l}])
      ok=false;
   end
end
if ~ok, error('Too many dual technique observations defined, aborting...');, end

%% Simulate the vertical and horizontal velocity field 
%
% We use the Matlab function peaks to simulate the vertical velocity field. 
% The horizontal velocity field is computed from the gradients. Peaks is
% defined over the xy range [-3 3], but we compute it over a slightly
% larger range. The grid points are defined by the vector u.

%u=-3.5:.1:4.5;
u=-3:.1:3;

% simulate vertical velocity field

[x,y,vz]=peaks(u);
[nx,ny]=size(vz);

% simulate height field 

hgt = 4+(x+y)/8+(x/8).^2+y.*exp(-(x/3).^2- y.^2);

% To obtain realistic x and y coordinates a scale factor and offset 
% is applied. Units of coordinates are in km.

scale=3;
offsetxy=[4 14]*scale; 

x=offsetxy(1)+scale*x;
y=offsetxy(2)+scale*y;

% horizontal velocity field

[vx,vy] = gradient(-1*vz,.2,.2);

%% make 3-D surf plot

figure('Name','Surface plots','Position',[680 558 840 420])

subplot(1,2,1)
surf(x,y,hgt)
xlabel('X [km]')
ylabel('Y [km]')
zlabel('hgt [m]')
title('Height')

subplot(1,2,2)
surf(x,y,vz)
xlabel('X [km]')
ylabel('Y [km]')
zlabel('vz [mm/y]')
title('Vertical velocity field')

%% make contour plot with quivers for the horizontal velocity

figure('Name','Contour plots','Position',[680 558 1260 420])

subplot(1,3,1)
contour(x,y,vz,'showtext','on')
title('Vertical and horizontal velocity field [mm/y]')
axis equal
hold on
quiver(x,y,vx,vy)
xlabel('X [km]')
ylabel('Y [km]')

% contour plots of the horizontal velocity field

subplot(1,3,2)
contour(x,y,vx,'showtext','on')
axis equal
xlabel('X [km]')
ylabel('Y [km]')
title('Horizontal velocity field (x-component) [mm/y]')

subplot(1,3,3)
contour(x,y,vy,'showtext','on')
axis equal
xlabel('X [km]')
ylabel('Y [km]')
title('Horizontal velocity field (y-component) [mm/y]')

    
%% Simulate observation points for the different techniques
%
% Observation points are randomly selected from the simulated grid, and 
% then assigned to the various observation techniques.

% Compute the totals

numtechniques=numel(techniqueId);   % Number of techniques
numdualtech=sum(npdualtech(:));    % Number of points with dual observations

% np is the total number of points

np=sum(nptech)-npalltech*(numtechniques-1)-numdualtech;

% Draw points randomly from integer grid cells. Since there can be doubles
% we initially generate more than strictly needed. The vector ix and iy
% are indices to x, y, and vx, vy and vz.

ix=randi(nx,np*2,1);
iy=randi(ny,np*2,1);

% Remove doubles and trim to total number of points

ind = sub2ind([nx,ny],ix,iy);
ind=unique(ind,'stable');
[ix,iy] = ind2sub([nx,ny],ind(1:np));

% Fill arrays with point coordinates and velocities

pntNo =[1:np]';
pntId=cell(np,1);
for k=1:np
  pntId{k}=sprintf('p%03d%03d',ix(k),iy(k));
end
pntCrd = [ y(iy,1) x(1,ix)' hgt(sub2ind([nx,ny],ix,iy)) ];
pntVel = [ vy(sub2ind([nx,ny],ix,iy)) vx(sub2ind([nx,ny],ix,iy)) vz(sub2ind([nx,ny],ix,iy)) ];

% Convert point coordinates to latitude and longitude, and change order in North, East and Up

pntYx=pntCrd(:,1:2)-repmat(mean(pntCrd(:,1:2)),[np 1]);
pntCrd=neu2plhsp([pntYx.*1000 pntCrd(:,3)],[53 6 0]);

% Generate "observation" matrix Ap of size np x numtechniques with zeros and ones 

Ap=zeros(np,numtechniques);
ip=0;
for k=1:npalltech
    ip=ip+1;
    Ap(ip,:)=ones(1,numtechniques);
end
for l2=1:size(npdualtech,2)
    for l1=l2+1:size(npdualtech,1)
       for k=1:npdualtech(l1,l2)
          ip=ip+1; 
          Ap(ip,[l2 l1])=[ 1 1 ];
       end
    end
end
for l=1:numtechniques
   for k=1:nptech(l)-sum(npdualtech(:,l))-sum(npdualtech(l,:))-npalltech
      ip=ip+1; 
      Ap(ip,l)=1;
   end
end

if ( ip ~= np ), error('Hey, this should not happen, ip should be equal to np.'); end

%% Plot the simulated field with points

figure('Name','Simulated point field')

contour(u,u,vz,'showtext','on')
hold on
quiver(u,u,vx,vy)
for k=1:numtechniques
  idxp=find(Ap(:,k));
  h(k)=plot(u(ix(idxp)),u(iy(idxp)),plotSymbol{k},'MarkerSize',7,'MarkerFaceColor','green');
end
legend(h,datasetId)
title('Vertical and horizontal velocity field (with points)')

%% Simulate observation epochs for the different techniques

% Make array with the epoch Id (and delta t)

nepochs=numel(t);
epochId=cellstr(num2str(t'))';

dt=0:nepochs-1;

% Generate "observation" matrix At of size numtechniques x nepochs  with zeros and ones 
%
% - simulate a delay (at most maxStartDelay) at the start of datasets
% - simulate a early end (at most maxEarlyEnd) at the end of datasets

ifirst=1+randi([0 maxStartDelay],numtechniques,1);
ilast=nepochs-randi([0 maxEarlyEnd],numtechniques,1);
if ~any(ifirst == 1), ifirst(randi(numtechniques,1,1))=1; end
if ~any(ilast == nepochs), ilast(randi(numtechniques,1,1))=nepochs; end

At=zeros(numtechniques,nepochs);
for l=1:numtechniques
  At(l,ifirst(l):ilast(l))=ones(1,ilast(l)-ifirst(l)+1);
end

%% Generate datasets
%
% In this section we first define the meaurement epochs, and then generate
% for each technique the space time observation matrix Y. The space time
% matrices are stored in the cell array obs. 

% Meta data for the datasets

creationDate=datestr(now,30);
projectFile=[ simid '_filenames.txt' ];
fileFormat='mat';

% Create structure for each dataset with space time matrix 

datasets={};
for k=1:numtechniques
  
  idxp=find(Ap(:,k));
  idxt=find(At(k,:));
  npi=numel(idxp);
  nti=numel(idxt);

  % Create stm structure 

  st = stm(simid,datasetId{k},techniqueId{k});
  % st = stm(simid,datasetId{k},techniqueId{k},'OBJECT');

  % Point and epoch information 

  st.numPoints=npi;
  st.numEpochs=nti;
  st.pntName=pntId(idxp);
  st.pntCrd=pntCrd(idxp,:);
  st.epochDyear=t(idxt);

  % Add harmonized pntId and epochId, and set projectFile name and
  % creationDate(projectId is set upon initialization).

  % subindexing is not possible with matfile objects. The code below
  % works for stuctures, but not for objects
  %
  % st.pntAttrib.pntId=pntId(idx);
  % st.epochAttrib.epochId=epochId;
  % st.datasetAttrib.projectFile=projectFile;
  % st.datasetAttrib.projectFileDate=creationDate;
  %
  % For objects use the following code instead
  
  pntAttrib=[];
  epochAttrib=[];
  pntAttrib.pntId=pntId(idxp);
  epochAttrib.epochId=epochId(idxt);

  datasetAttrib=st.datasetAttrib;
  datasetAttrib.projectFile=projectFile;
  datasetAttrib.projectFileDate=creationDate;
  st.datasetAttrib=datasetAttrib;

  % Add the technique dependent data
  
  switch techniqueId{k}
    case 'lev'     
       % Leveling observations with arbitrary height reference point
       % (simulating a fixed reference benchmark)
       Y = pntVel(idxp,3)*dt(idxt);
       Y = Y - repmat(Y(1,:),[size(Y,1) 1]);
       st.obsTypes = { 'Up' };
       st.obsData = Y;
       st.sensitivityMatrix = [ zeros(npi,2) ones(npi,1) ];
       % create stochastic model (block diagonal)
       C=stmstochmodel(stochModelLev{1},[],pntYx(idxp,1:2),t(1));
       stochData=zeros(npi,npi,nti);
       for l=1:nti
          stochData(1:npi,1:npi,l)=C;
       end
       st.stochModel{1}={['covmatrix(format=blkdiag,nd=' num2str(nti) ')']};
       st.stochData=stochData;
    case 'gnss'
       Y=nan(npi,nti,3);
       Y(:,:,1)= pntVel(idxp,1)*dt(idxt);
       Y(:,:,2)= pntVel(idxp,2)*dt(idxt);
       Y(:,:,3)= pntVel(idxp,3)*dt(idxt);
       st.obsTypes = { 'North' 'East' 'Up' };
       st.obsData = Y;
       st.sensitivityMatrix=nan(npi,3,3);
       st.sensitivityMatrix(:,:,1) = [ ones(npi,1) zeros(npi,1) zeros(npi,1) ];
       st.sensitivityMatrix(:,:,2) = [ zeros(npi,1) ones(npi,1) zeros(npi,1) ];
       st.sensitivityMatrix(:,:,3) = [ zeros(npi,1) zeros(npi,1) ones(npi,1) ];     
       st.stochModel=stochModelGNSS;
       st.stochData=[];
     case 'insar'
       % Double difference SAR observations
       Y = -1.0 * ( sind(incAngle(k)) * ( ...
                       cosd(azAngle(k))*pntVel(idxp,1) +      ...
                       sind(azAngle(k))*pntVel(idxp,2)   )  + ...
                    cosd(incAngle(k))*pntVel(idxp,3)   ...
                  ) * dt(idxt);
       Y = Y - repmat(Y(:,1),[1 size(Y,2)]);
       Y = Y - repmat(Y(1,:),[size(Y,1) 1]);
       st.obsTypes = { 'los' };
       st.obsData = Y;
       st.sensitivityMatrix = -1.0 * ones(npi,1) * [ cosd(azAngle(k))*sind(incAngle(k)) sind(azAngle(k))*sind(incAngle(k)) cosd(incAngle(k)) ];
       % cannot use subindexing for objects
       %st.pntAttrib.incAngle = ones(npi,1)*incAngle(k);
       %st.pntAttrib.azAngle = ones(npi,1)*azAngle(k);
       pntAttrib.incAngle = ones(npi,1)*incAngle(k);
       pntAttrib.azAngle = ones(npi,1)*azAngle(k);
       st.stochModel=stochModelInSAR;
       st.stochData=[];
     otherwise
       error('unknown technique')
  end
  % save attributes at the last (objects don't support subindexing)
  st.epochAttrib=epochAttrib;
  st.pntAttrib=pntAttrib;
  datasets{k}=st;
end


%% Randomly remove some observations to simulate missing observations
%
% * for levelling remove some complete epochs and random observations
% * for GPS, remove for some points, starting epochs
% * for InSAR, keep all (for the time being)

if removeobs

for k=1:numtechniques

  npi=datasets{k}.numPoints;
  nti=datasets{k}.numEpochs;
 
  switch techniqueId{k}
    case 'lev'
       % remove random observables
       % - percMissingLev is the percentage of points with missing observations
       % - maxMissingLev is the maximum number of missing observations for
       %   points with missing observations

       numMissing=randi([1 maxMissingLev],npi,1);
       numMissing(randperm(npi,round(npi*(100-percMissingLev)./100)))=0;    
%      sum(numMissing>0)
%      mean(numMissing)
%      max(numMissing)
       
       Y=datasets{k}.obsData;
       for l=1:npi
          if numMissing(l)==0, continue; end
          Y(l,randperm(nti,numMissing(l)))=nan;
       end
       
       % remove numMissingEpochsLev epoch(s)      
       
       iremove=randperm(nti,numMissingEpochsLev);    
       
       datasets{k}.numEpochs=nti-numMissingEpochsLev;
       datasets{k}.epochDyear(iremove)=[];
       epochAttrib=datasets{k}.epochAttrib;
       epochAttrib.epochId(iremove)=[];
       datasets{k}.epochAttrib=epochAttrib;
       
       Y(:,iremove)=[];
       datasets{k}.obsData=Y;       

       datasets{k}.stochModel={['covmatrix(format=blkdiag,nd=' num2str(nti-numMissingEpochsLev) ')']};
       datasets{k}.stochData(:,:,iremove)=[];
       
    case 'gnss'
       % remove observations from the start
       % - percMissingGNSS is the percentage of points with missing observations
       % - maxMissingGNSS is the maximum number of missing observations for
       %   points with missing observations
     
       numMissing=randi([1 maxMissingGNSS],npi,1);
       numMissing(randperm(npi,round(npi*(100-percMissingGNSS)./100)))=0;    

       Y=datasets{k}.obsData;
       for l=1:npi
          if numMissing(l)==0, continue; end
          Y(l,1:numMissing(l),:)=nan;
       end
       datasets{k}.obsData=Y;

    case 'insar'
    otherwise
       error('unknown technique')
  end
end

end

%% Plot the 'truth' observations

figure('Name','Truth observations','Position',[50 50 1360 860])

kk=1;
for k=1:numel(datasets)
   obsTypes=datasets{k}.obsTypes;
   for l=1:numel(obsTypes)
      %figure
      subplot(3,3,kk);kk=kk+1;
      Y=datasets{k}.obsData(:,:,l);
      %imagesc(Y);
      %cmax=nanmax(nanmax(abs(Y)));
      %caxis([-cmax cmax])
      %colorbar
      waterfall(Y)
      title([ datasets{k}.datasetId ' (' obsTypes{l} ')'] )
   end
end

%% Save the 'truth' observation datasets

fid=fopen([ simid '_filenames_truth.txt' ],'w');
for k=1:numel(datasets)
   datasetFilename=stmwrite(datasets{k},[ simid '_' datasetId{k} '_truth.mat']);
   fprintf(fid,'%s\n',datasetFilename);
end
fclose(fid);


%% Simulate observation noise

for k=1:numtechniques
  
  st=datasets{k};
  idxp=find(Ap(:,k));
  %idxt=find(At(k,:));
  
  Y=st.obsData;
  npi=size(Y,1);
  nti=size(Y,2);

  numObsTypes=numel(st.obsTypes);
%  if numObsTypes == 1
%    R=chol(stmstochmodel(st.stochModel{1},st.stochData,pntYx(idx,1:2),st.epochDyear));
%    Y(:) = Y(:) + R * randn(npi*nti,1);  
%  else
    for l=1:numObsTypes
      if strcmpi(stochModel,'simple')
         Y(:,:,l) = Y(:,:,l) + randn(npi,nti);
      else
         R=chol(stmstochmodel(st.stochModel{l},st.stochData,pntYx(idxp,1:2),st.epochDyear));
         Y(:,:,l) = Y(:,:,l) + reshape( R * randn(npi*nti,1), [npi,nti]);
      end
    end
%  end
  datasets{k}.obsData=Y;

end

%% Plot the 'noisy' observations

figure('Name','Noisy observations','Position',[50 50 1360 860])

kk=1;
for k=1:numel(datasets)
   obsTypes=datasets{k}.obsTypes;
   for l=1:numel(obsTypes)
      %figure
      subplot(3,3,kk);kk=kk+1;
      Y=datasets{k}.obsData(:,:,l);
      %imagesc(Y);
      %cmax=nanmax(nanmax(abs(Y)));
      %caxis([-cmax cmax])
      %colorbar
      waterfall(Y)
      title([ datasets{k}.datasetId ' (' obsTypes{l} ')'] )
   end
end

%% Save the observation datasets (with noise included)

fid=fopen([ simid '_filenames.txt' ],'w');
for k=1:numel(datasets)
   datasetFilename=stmwrite(datasets{k});
   fprintf(fid,'%s\n',datasetFilename);
end
fclose(fid);

%% Generate dataset with true displacement field

% Create stm structure 

st = stm(simid,'truth','displ');

% Point and epoch information 

st.numPoints=np;
st.numEpochs=nepochs;
st.pntName=pntId;
st.pntCrd=pntCrd;
st.epochDyear=t;

% Add harmonized pntId and epochId, and set projectFile name and
% creationDate(projectId is set upon initialization).
  
pntAttrib=[];
epochAttrib=[];
pntAttrib.pntId=pntId;
epochAttrib.epochId=epochId;
st.pntAttrib=pntAttrib;
st.epochAttrib=epochAttrib;

datasetAttrib=st.datasetAttrib;
datasetAttrib.projectFile=projectFile;
datasetAttrib.projectFileDate=creationDate;
st.datasetAttrib=datasetAttrib;

% True displacement
  
Y=nan(np,nepochs,3);
Y(:,:,1)= pntVel(:,1)*dt;
Y(:,:,2)= pntVel(:,2)*dt;
Y(:,:,3)= pntVel(:,3)*dt;
st.obsTypes = { 'North' 'East' 'Up' };
st.obsData = Y;
st.sensitivityMatrix=nan(np,3,3);
st.sensitivityMatrix(:,:,1) = [ ones(np,1) zeros(np,1) zeros(np,1) ];
st.sensitivityMatrix(:,:,2) = [ zeros(np,1) ones(np,1) zeros(np,1) ];
st.sensitivityMatrix(:,:,3) = [ zeros(np,1) zeros(np,1) ones(np,1) ];     
st.stochModel{1}={ 'WN(sigma=0)','spatialcov(rho=Inf)' };  
st.stochData=[];

% Save the truth dataset

stmwrite(st);

% Finish the function

fprintf('%s finished at %s  (elapsed time %.2f s)\n',progname,datestr(now),toc);

end
% [End of Function]