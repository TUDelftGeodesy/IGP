function st=stmread(stmfile,mode,varargin)
%READSTM    Read Space Time Matrix dataset from file.
%   ST=STMREAD(STMFILE) reads the Space Time Matrix structure ST from
%   the the file STMFILE. 
%
%   ST=STMREAD(STMFILE,'NODATA') reads the Space Time Matrix structure
%   ST but defering the actual read of the observation, auxiliary and 
%   stochastic data. This option is useful for large files that would
%   otherwise not fit into internal memory. Use another call to STMREAD
%   to read (a subset of) the observations, auxiliary or stochastic data
%   into memory.
%
%   ST=STMREAD(STMFILE,'NOSTOCHDATA') idem deferring only stochastic data
%   (useful for full co-variance matrices).
%
%   ST=STMREAD(STMFILE,'OBJECT') reads the Space Time Matrix as a dynamic
%   object using Matlab matfile utility. This is only supported for Matlab
%   V7.3 matfiles and higher.
%
%   MAT=STMREAD(STMFILE,DATATYPE), with DATATYPE either 'OBSDATA', 'AUXDATA'
%   or 'STOCHDATA', reads the observation, auxiliary or stochastic data 
%   from STMFILE. The data is returned as a Matlab matrix MAT (Space Time
%   Matrix).
%
%   MAT=STMREAD(STMFILE,DATATYPE,PNTIDX,EPOIDX,COMPIDX) reads a subset
%   of DATATYPE. The result is returned as a Matlab matrix MAT (space Time
%   Matrix) for the requested subset. The range of the indices is not
%   checked in this function. If indices are used, PNTIND and EPOIDX both
%   must be specified explicitly. If, COMPIDX is missing, a value of 1 (first
%   component) is assumed. This option is only supported for NetCDF files 
%   (not yet implemented) and mat V7.3 files.
%
%   MAT=STMREAD(STMOBJECT,DATATYPE,PNTIDX,EPOIDX,COMPIDX) reads a subset
%   of DATATYPE from a matfile object STMOBJECT, created by an earlier
%   call to STMREAD.
%
%   The STM file is either stored as a Matlab V7.3 mat file or as a 
%   Cupido V2 NetCDF file. Support for NetCDF has not yet been implemented.
%
%   The Space Time Matrix structure ST has the following fields
%
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
%      .stochModel      Stochastic model description    cell array of strings
%      .stochData       Stochastic model data           single [numPoints,numEpochs,dimStochData]
%      .inputDatasets() Structure array with for each input dataset  struct array 
%
%   For creating Space Time Matrix datasets and more details on the dataset
%   stucture see the help of the STM function.
%
%   Examples:
%
%      st=stmread('simtrue1_sarAsc1.mat')
%      st=stmread('simtrue1_sarAsc1.mat','nodata')
%      obs=stmread('simtrue1_sarAsc1.mat','obsdata')
%      obs=stmread('simtrue1_sarAsc1.mat','obsdata',1:10,1:5)
%      obs=stmread('simtrue1_sarAsc1.mat','obsdata',1:10,1:5)
%      obs=stmread('simtrue1_gps.mat','obsdata',1:5,3:12,1:3)
% 
%      st=stmread('simtrue1_sarAsc1.mat','object')
%      obs=stmread(st,'obsdata')           % same as st.obsdata
%      obs=stmread(st,'obsdata',1:10,1:5)  % same as st.obsdata(1:10,1:5)
%
%   See also STM and STMWRITE.
%
%   (c) Delft University of Technology, 2020.

% Created:  20 August 2020 by Hans van der Marel
% Modified: 25 August 2020 by Hans van der Marel
%           - Added support for reading subsets
%           - Added examples to the help
%           26 Oct 2023 by Hans van der Marel
%           - Added NOSTOCHDATA option

% Check input arguments

if nargin <1, error('Need to have at least the name of the space time matrix file'); end
if ischar(stmfile) || isstring(stmfile)
   if ~exist(stmfile,'file')
      error(['Input space time matrix file ' stmfile ' does not exist.'])
   end
   [~,~,ext]=fileparts(stmfile);
elseif isobject(stmfile) 
   ext='.mat';
else
   error('First input argument must be a character string or object.');
end

if nargin < 2
  mode='DATA';
end

dosubset=false;
if nargin > 2 && nargin <= 5 && any(strcmpi({'OBSDATA' 'AUXDATA' 'STOCHDATA'},mode))
  if nargin < 4 
     error('PNTIDX and EPOIDX  must both be supplied')
  elseif nargin == 5
     compidx=varargin{3};
  else
     compidx=1;
  end
  pntidx=varargin{1};
  epochidx=varargin{2};
  dosubset=true;
elseif nargin > 2 && nargin <= 5
  warning('PNTIDX, EPOIDX and COMPIDX only needed when MODE is a DATATYPE.') 
elseif nargin > 5
  error('Too many input arguments.')  
end

% Read the space time matrix file


switch lower(ext)
    case '.mat'
        % Read space time matrix from Matlab mat file
        switch upper(mode)
            case { 'DATA' }
               st=load(stmfile);
            case { 'NODATA' }
               % Only load variables that do not match "Data" at end of name
               st=load(stmfile,'-regexp','...(?<!Data)$');
            case { 'NOSTOCHDATA' }
               % Only load variables that do not match "StochData" at end of name
               st=load(stmfile,'-regexp','...(?<!stochData)$');
            case { 'OBJECT' 'DYNAMIC' }
               st=matfile(stmfile);
            case { 'OBSDATA' 'AUXDATA' 'STOCHDATA' }
               var=strrep(lower(mode),'data','Data');
               if dosubset
                  if isobject(stmfile)
                      st=stmfile.(var)(pntidx,epochidx,compidx);
                  else
                      % create temporary object for stmfile
                      stmobj = matfile(stmfile);
                      st=stmobj.(var)(pntidx,epochidx,compidx);
                  end
               else
                  st=load(stmfile,var);
                  st=st.(var);
               end
            otherwise
               error(['Wrong 2nd argument ' mode ' (must be DATA, NODATA, NOSTOCHDATA, OBSDATA, AUXDATA, STOCHDATA or OBJECT'])
        end
    otherwise
        % Read space time matrix from NetCDF file
        % stm=stmreadnetcdf(stmfile,mode)  % Internal function 
        error('NetCDF file format not yet supported')
end

% Add filename to the output structure (this fails for objects)

%st.currentFileName=stmfile;

end