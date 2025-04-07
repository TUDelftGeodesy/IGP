function st=stm(projectId,datasetId,techniqueId,mode,filename)
%STM    Create a Space Time Matrix dataset structure or object.
%   ST=STM() creates an empty Space Time Matrix dataset structure ST.
%   Some of the dataset attributes (createdBy, creationDate, ...) in ST 
%   are set with values determined from the OS. 
%
%   ST=STM(PROJECTID,DATASETID,TECHNIQUEID) creates a Space Time Matrix 
%   dataset structure ST with PROJECTID, DATASETID and TECHNIQUEID.
%   The TECHNIQUEID is checked against allowable values and is used
%   to set already some default values.
%
%   ST=STM(...,'OBJECT') creates a writable Space Time Matrix dataset 
%   matfile object PROJECTID_DATASETID.mat in the current directory.
%   With ST as matfile object variables can be accessed and changed 
%   without loading the file into memory. You can load or save parts 
%   of variables. Partial loading and saving of variables requires 
%   less memory than the load and save commands. Access to variables 
%   is similar to accessing fields of structs. 
%
%   ST=STM(...,'OBJECT',FILENAME) creates a writable Space Time Matrix 
%   dataset matfile object FILENAME. FILENAME can include a full or 
%   partial path, otherwise matfile searches along the MATLAB path.
%   If FILENAME is not given a default filename is created from
%   the PROJECTID and DATASETID.
%
%   ST=STM('OBJECT',FILENAME) creates a writable Space Time Matrix 
%   dataset matfile object FILENAME. FILENAME can include a full or 
%
%   'OBJECT' one of the possible values of MODE. Allowable values
%   for  MODE are
%
%       'STRUCT', 'DATA' : ST is a structure. The structure can be
%             saved as mat-file or NetCDF file when the structure
%             is completed using stmwrite. 
%       'OBJECT', 'DYNAMIC' : ST is a matfile object. Object variables 
%             can be accessed and changed without loading the file 
%             into memory. Access to variables is similar to accessing 
%             fields of structs. 
%       'NODATA' : ST is a structure. Does the same as 'STRUCT', but
%             provided for compability with stmread. 
%
%   Future MODEs include dynamic access to NetCDF files. The default
%   for MODE is 'STRUCT'.  
%
%   See also STMREAD, STMWRITE, STMSTOCHMODEL and MATFILE.
%
%   The STM file is either stored as a Matlab V7.3 mat file or as a 
%   Cupido V2 NetCDF file.
%
%   The Space Time Matrix structure has the following fields
%
%   FieldName       Description                     Type                       Values/Examples/Notes
%   --------------  ------------------------------- ------------------------   -----------------------
%   datasetId       Solution Identifier             char string (free)         { 'lev1'  'gps1'  'sarDsc1'  'sarAsc1' 'sarAsc2'} 
%   techniqueId     Technique Identifier            char string (reserved)     { 'lev'  'gnss'  'insar' 'displ' }
%
%   datasetAttrib   Dataset Attributes              struct                     (see note 9)
%      .createdBy      Username of the account running the software
%      .creationDate   File creation date in ISO format
%      .createdAs      File name at the time the file was created, with full path
%      .softwareName   Name of the software that created the dataset/file
%      .softwareOptionsOptions that were used to create the dataset/file
%      .fileFormat     File format
%      .fileName       File name at the time the file was read
%      .projectId      Project Id or name                                      (see note 10) 
%      .projectFile    Project file name                                       (see note 2,10)
%      .projectFileDate  Creation/modification date of the project file        (see note 2,10)   
%                        in ISO format
%   inputDatasets()  Structure array with for each input dataset  struct array (see note 11)
%      .datasetId
%      .techniqueId
%      .datasetAttrib
%      .numPoints
%      .numEpochs
%      .inputDatasets
%   techniqueAttrib Technique Attributes            struct 
%      .system         System name                     char string             { 'GPS' 'GNSS' 'TSX 'Sentinel-1' 'Radarsat-2'  ... }
%      .mode           Technique specific mode         char string             { 'CORS' 'Campaign'  'Asc'  'Dsc' }
%      ...             Other attributes (free names)   char string 
%                      e.g. frequency, track, ...
%   globalAttrib    Global Attributes               struct                     (see note 12)
%                      NetCDF globalattributes
%
%   numPoints       Number of Points                int scalar
%   numEpochs       Number of Epochs                int scalar
%
%   pntName         Point name                      cell array [numPoints]     (see note 1) 
%   pntCrd          Point Cooordinates (lat,lon,h)  double [numPoints,3] matrix
%   pntAttrib       Point Attributes ....           table  [numPoints,* ]      (see note 3, 13)
%      .pntId          Harmonized Point Identifier     cell array [numPoints]  (see note 1,2)               
%      .pntClass       Point class (conform CupidoV1)  cell array [numPoints]  (see note 4)
%      .incAngle       Incidence angle [deg] (InSAR)   single [numPoints,1] matrix
%      .azAngle        Azimuth angle [deg] (InSAR)     single [numPoints,1] matrix
%
%   epochDyear    Decimal year                      double [numEpochs] matrix 
%   epochAttrib   Epoch Attributes                  table [numEpochs,*]        (see note 3, 13)
%      .epochId      Harmonized Epoch identifier       cell array [numEpochs]  (see note 1,2) 
%      .prjName      Project name (conform CupidoV1)   cell array [numEpochs]  (see note 4)
%      .prjClass     Project class (conform CupidoV1)  cell array [numEpochs]  (see note 4)
%
%   obsTypes      Observation Types                 cell array  (reserved)     { 'East' 'North' 'Up' 'los' }        (see note 5)
%   parTypes      Parameter Types                   cell array  (reserved)     { 'North' 'East' 'Up' 'Compaction' } (see note 5,6)
%   auxTypes      Auxiliary data Types              cell array (free)          { 'amp' 'pow' 'snr' ...  }           (see note 5)
%
%   obsData       Observation Space Time Matrix     single [numPoints,numEpochs,numObsTypes]   (see note 5,7)
%   auxData       Auxiliary data Space Time Matrix  single [numPoints,numEpoch,numAuxTypes]    (see note 5,7)
%
%   sensitivityMatrix  Sensitivity Matrix           single [numPoints,numParTypes,numObsTypes]      (see note 5)
%
%   stochModel    Stochastic model description      cell array of char strings (see note 8)
%   stochData     Stochastic model data             single [numPoints,numEpochs,dimStochData] (see note 7)
%
%
%   Note 1: pntId and epochId are harmonized point and epoch identifiers,
%      meaning different datasets use the same identifier so that we can
%      easily match points and epochs. The coordinates pntCrd and epoch times 
%      epochDyear refer to the actual point and epoch (they are not harmonized),
%      as is pntName which contains the original dataset dependent point
%      names (or the harmonized variant). The harmonization requires
%      an extra processing step to prepare pntId and epochId, which are
%      matching (harmonized) character strings between the datasets. 
%      Therefore, pntId and/or epochId may, or may not, be present; hence
%      the are part of pntAttrib which contains optional items.
%
%   Note 2: If the pntId and/or epochId field is not present in the pntAttrib,
%      then, for instance the integration module, must use the coordinates 
%      and/or epoch times are used to identify matching points and epochs 
%      between datasets (using simple algorithms).
%
%   Note 3: All pntAttrib and epochAttrib fields are optional and not restricted 
%      to the examples given in the description. Also note that pntAttrib 
%      and epochAttrib are not attributes in NetCDF sense; they have to
%      be stored as numeric and/or string arrays in NetCDF.
%
%   Note 4: pntClass, prjName and prjClass from Cupido V1 files are
%      stored in (optional) pntAttrib and epochAttrib fields.
%
%   Note 5: The parameter and observation dimension is not given explicitly in the 
%      datastructure. They follow simply from the size of the obsTypes, parTypes 
%      and auxTypes arrays, with numXxxTypes=numel(xxxTypes),
%
%        numObsTypes   Number of observed components      int scalar     1,2,3 
%        numParTypes   Displacement dimension             int scalar     3 or 4
%        numAuxTypes   Number of layers in auxiliary STM  int scalar     1,2,3,...
%
%      We suggest to use these names in the software to improve readability.
%
%   Note 6: The order of the parameter components is fixed to North, East, Up !!!!
%      The dimension dimPar is therefore 3 or 4 (4 when the optional compaction
%      component is include), even when we estimate for instance only one 
%      component (e.g. height). The name for the fourth component, compaction, 
%      has to be decided. This will be only added when information on compaction 
%      is available from the observations.
%
%   Note 7: Reading of obsData, auxData and stochData can be deferred to
%      later and may include subsets of the data only. This is done to
%      - keep the size of the stm structure small (and in memory)
%      - allow for sequential processing of the observation data
%      Smaller datasets may be read into memory (default). For larger
%      datasets we may use matread for V7.3 mat files and/or a dedicated
%      stmreadobs function (NetCDF and V7.3 mat files)
%
%   Note 8: The stochastic model can be specified in terms of a) the co-variance
%      matrix, b) precision matrix (inverse of the covariance matrix) or c) the
%      square root (Choleski factor) of either of them. The matrix itself can
%      either be storead as full, lower triangle, block diagonal or as 
%      diagonals, or specified in a function format in stochModel. See
%      stmstochmodel function for details.
%
%   Note 9: datasetAttrib is a structure with dataset attributes. All fields,
%      except fileName and fileFormat, are derived from the file contents. 
%      However, fileName (with full path) is the name of the file at the
%      time of reading and comes from the o/s, and the fileFormat is known
%      to the software reading the file.
%
%   Note 10: pntId and epochId contain harmonized names between datasets.
%      This is for the integration module that needs to know the common
%      common points and epochs, and because there is a dedicated module
%      to prepare this information. This implies that we must have a 
%      method to establish that the pntId and epochId used in different
%      datasets indeed use corresponding identifiers (and not that they are
%      from different projects). For this we require that the projectId,
%      projectFile and projectFileDate fields are identical for different
%      datasets. The projectFile is the name with the full path of the file 
%      that is created by the module to identify common points and epoch,
%      and has been used to generate the pntId and epochId. The content
%      of this file is not important: only the name and filedate (ISO format)
%      must be the same.
%      If the projectFile name is empty, or not given, then this is a
%      signal that the points and epochs have not been harmonized.
%
%   Note 11: inputDatasets is a structure array with information on the
%      input datasets that have been used to create this (output) dataset.
%      It contains the datasetId, techniqueId, datasetAttrib, numPoint,
%      numEpochs and inputDatasets from the input datasets. Including
%      inputDatasets from the input datasets makes a full traceback
%      possible of all the previous computation steps. If the structure
%      array the traceback ends. 
%
%   Note 12: NetCDF global attributes are stored in globalAttrib as
%      a structure. For NetCDF this has to be converted from/to a cell
%      array. The Global attributes (c.f. NETCDF CF-1.6) are defined by a cell
%      array with two columns. Mandatory attributes are
%
%      globalattributes = { ...
%        'title'        , 'A succinct description of what is in the dataset.'  ; ...
%        'institution'  , 'Specifies where the original data was produced.' ; ...
%        'source'       , 'The method of production of the original data. If it was model-generated, source should name the model and its version, as specifically as could be useful. If it is observational, source should characterize it, e.g., surface observation or radiosonde.' ; ...
%        'history'      , 'Provides an audit trail for modifications to the original data. Well-behaved generic netCDF filters will automatically append their name and the parameters with which they were invoked to the global history attribute of an input netCDF file. We recommend that each line begin with a timestamp indicating the date and time of day that the program was executed.' ; ...
%        'references'   , 'Published or web-based references that describe the data or methods used to produce it.' ; ...
%        'comment'      , 'Miscellaneous information about the data or methods used to produce it.' ; ...
%        'Conventions'  , 'CF-1.6' ; ...
%         'email'        , 'h.vandermarel@tudelft.nl' ; ...
%        'version'      , ' ' ; ...
%        'terms_for_use', 'These data can be used freely for research purposes provided that the following source is acknowledged: Actueel Hoogtebestand Nederland (AHN): Rijkswaterstaat + WaterSchapsHuis' ; ...
%        'disclaimer'   , 'This data is made available in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.' ; ...
%        'featureType'  , 'Specifies the type of discrete sampling geometry to which the data in the file belongs, and implies that all data variables in the file contain collections of features of that type.' ; ... 
%      };
%
%      The feature type is 'timeSeries' , this is the standard for all datasets 
%
%   Note 13: pntAttrib and epochAttrib may either be a table or a structure.
%      Table is supported from Matlab 2006 upwards and are convenient 
%      because it guaranties all arrays have the same length. A structure
%      is a bit more flexible and supported by older Matlab versions.
%      The Matlab syntax is more or less the same for table and structures.
%
%   (c) Delft University of Technology, 2020.

% Created by Hans van der Marel
%
% Version history:
%
%    v0.1   15 July 2019    Initial release IspSimDefo2.m
%    v0.2   16 April 2020   Major update
%    v0.3   10 August 2020  Improve consistency and naming, stochastic model added 
%    v0.3-1 13 August 2020  Further iteration based on discussion with Freek 
%    v0.3-2 19 August 2020  Split space time matrix in geometric and non
%                           geometric components; geometric components are
%                           alligned with the sensitivity matrix
%    v0.4   20 August 2020  Birth of stmread with many updates
%    v0.5   14 Sept 2020    Updates for stochastic model
%    v0.6   25 Sept 2020    Small change to avoid . indexing on functions
%                           (to enable use of older Matlab versions) 
%    v0.6-1 01 Oct 2020     Added 'displ' as new "technique" for displacements
%                           (to enable use of older Matlab versions) 
%    v0.6-2 20 Oct 2021     Added softwareOptions to datasetAttrib and
%                           solved bug for the stm() call
%
% Open issues:
%    store input options in data structure

%   Mapping of Cupido(v1) variables to Cupido-STM format variables
%
%      Cupido(v1)       stm structure format 
%      ---------------- --------------------------------------------------
%      globalattributes globalAttrib 
%
%      pntname          pntName
%      pntcrd           pntCrd 
%      pntclass         pntAttrib.class 
%
%      prjname          epochAttrib.prjName
%      prjepoch         epochDyear            [Matlab date number -> Dyear]
%      prjclass         epochAttrib.class
%
%      obstable         observation table with index to from_point, to_point 
%                       and project [only possible if from_point is the
%                       same for a single epoch]
%      sdobs            vector with the observed height difference [m]
%      sdcov            covariance matrix [m]
%      sdobsflag        integer observation flag (default 0)
%      sensitivity      sensitivity matrix [0-1] (default [ 0 0 1 ])
%      epoch            array with epoch [Matlab date number] (default prjepoch)

 
% Check input arguments and provide default values 

if nargin == 1 || nargin > 5
    error('this function expects 0, 2, 3, 4 or 5 input arguments')
end
if nargin == 0 || nargin == 3 
    mode='STRUCT';
end
if nargin == 2
    mode=projectId;
	filename=datasetId;
end
if nargin <= 2
    projectId='';
    datasetId='';
    techniqueId='';
end

fileFormat='mat';
if nargin ~= 2 && nargin ~= 5
    filename=[projectId '_' datasetId '.' fileFormat ];
end

% Check value for techniqueId

techniqueId=lower(techniqueId);
switch techniqueId
    case 'lev'
        obsTypes={'Up'};
        parTypes={'North' 'East' 'Up'};
    case {'gnss','displ'}
        obsTypes={'North' 'East' 'Up'};
        parTypes={'North' 'East' 'Up'};
    case 'insar'
        obsTypes={'los'};
        parTypes={'North' 'East' 'Up'};
    otherwise
        if ~isempty(techniqueId)
            error('invalid techniqueId')
        end
        obsTypes=[];
        parTypes=[];
end

% Create default datasetAttrib structure
%
%   datasetAttrib   Dataset Attributes              struct
%      .createdBy      Username of the account running the software
%      .creationDate   File creation date in ISO format
%      .createdAs      File name at the time the file was created, with full path
%      .softwareName   Name of the software that created the dataset/file
%      .fileFormat     File format
%      .fileName       File name at the time the file was read
%      .projectId      Project Id or name  
%      .projectFile    Project file name 
%      .projectFileDate  Creation/modification date of the project file    
%                        in ISO format

%createdBy=getenv('username');
createdBy={getenv('USER'),getenv('username'),'unknown'}; 
createdBy=createdBy{~cellfun('isempty',createdBy)};
creationDate=datestr(now,30);
dbstackCall=dbstack(1);
if ~isempty(dbstackCall)
   softwareName=dbstackCall.name;
else
   softwareName='';
end

datasetAttrib = struct( ...
    'createdBy',createdBy,...
    'creationDate',creationDate,...
    'createdAs',filename,...
    'softwareName',softwareName,...
    'softwareOptions',[],...
    'fileFormat',fileFormat,...
    'projectId',projectId,...
    'projectFile','',...
    'projectFileDate','' );

% Create empty inputDataset stucture array
% 
%   inputDatasets() Structure array with for each input dataset  struct array 
%      .datasetId
%      .techniqueId
%      .datasetAttrib
%      .numPoints
%      .numEpochs
%      .inputDatasets

inputDatasets=struct( ...
    'datasetId',{},...
    'techniqueId',{},...
    'datasetAttrib',{},...
    'numPoints',{},...
    'numEpochs',{},...
    'inputDatasets',{});

% Create space time matrix structure / object
%
%   st
%      .datasetId       Solution Identifier             char string (free)
%      .techniqueId     Technique Identifier            char string (reserved)
%      .datasetAttrib   Dataset Attributes              struct
%      .techniqueAttrib Technique Attributes            struct 
%      .globalAttrib    Global Attributes               struct
%      .numPoints       Number of Points                int scalar
%      .numEpochs       Number of Epochs                int scalar
%      .pntName         Point name                      cell array [numPoints] 
%      .pntCrd          Point Cooordinates (lat,lon,h)  double [numPoints,3] matrix
%      .pntAttrib       Point Attributes ....           table  [numPoints,* ] 
%      .epochDyear      Decimal year                    double [numEpochs] matrix 
%      .epochAttrib     Epoch Attributes                table [numEpochs,*] 
%      .parTypes        Parameter Types                 cell array  (reserved)
%      .obsTypes        Observation Types               cell array  (reserved)
%      .auxTypes        Auxiliary data Types            cell array (free)
%      .obsData         Observation Space Time Matrix   single [numPoints,numEpochs,dimObs]
%      .auxData         Auxiliary data Space Time Matrix  single [numPoints,numEpoch,dimAux]
%      .sensitivityMatrix  Sensitivity Matrix           single [numPoints,dimPar,dimObs]
%      .stochModel      Stochastic model description    cell array
%      .stochData       Stochastic model data           single [numPoints,numEpochs,dimStochData]
%      .inputDatasets() Structure array with for each input dataset  struct array 


% Optionally create dynamic object

switch upper(mode)
    case {'STRUCT' 'DATA' 'NODATA'}
        % Create a structure
        st = struct(  ...
            'datasetId',datasetId, ...
            'techniqueId',techniqueId,...
            'datasetAttrib',{ datasetAttrib },...
            'techniqueAttrib',[],...
            'numPoints',[],...
            'numEpochs',[],...
            'pntName',[],...
            'pntCrd',[],...
            'pntAttrib',[],...
            'epochDyear',[],...
            'epochAttrib',[],...
            'parTypes',{ parTypes },...
            'obsTypes',{ obsTypes },...
            'auxTypes',[],...
            'obsData',[],...
            'auxData',[],...
            'sensitivityMatrix',[],...
            'stochModel',[],...
            'stochData',[],...
            'inputDatasets',{ inputDatasets },...
            'globalAttrib',[]);
    case {'OBJECT' 'DYNAMIC' }
        % Create an matfile object (written instantly to disk) 
        st = matfile(filename,'Writable',true);
        st.datasetId=datasetId;
        st.techniqueId=techniqueId;
        st.datasetAttrib=datasetAttrib;
        st.techniqueAttrib=struct([]);
        st.numPoints=[];
        st.numEpochs=[];
        st.pntName=[];
        st.pntCrd=[];
        st.pntAttrib=struct([]);
        st.epochDyear=[];
        st.epochAttrib=struct([]);
        st.parTypes=parTypes;
        st.obsTypes=obsTypes;
        st.auxTypes=[];
        st.obsData=[];
        st.auxData=[];
        st.sensitivityMatrix=[];
        st.stochModel=[];
        st.stochData=[];
        st.inputDatasets=inputDatasets;
        st.globalAttrib=[];
    otherwise
        error('Invalid mode, use STRUCT / DATA / NODATA or OBJECT / DYNAMIC')        
end

end
