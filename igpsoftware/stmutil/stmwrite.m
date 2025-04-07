function varargout=stmwrite(st,stmfile,mode,varargin)
%STMWRITE    Write Space Time Matrix structure to file.
%   STMFILE=STMWRITE(ST) writes the Space Time Matrix structure ST 
%   to the file STMFILE that is set in the datasetAttrib of the structure.
%   The function returns the name of the STMFILE.
%
%   STMWRITE(ST,STMFILE) writes the Space Time Matrix structure ST 
%   to the file STMFILE, overruling the filename set inside the structure.
%
%   The STM file is either stored as a Matlab V7.3 mat file or as a Cupido 
%   V2 NetCDF file. Default is to use the extension from the STMFILE 
%   (.mat or .nc) to set the output format. Support for NetCDF is not 
%   yet implemented.
%
%   STMWRITE(ST,STMFILE,FMT) writes the Space Time Matrix structure 
%   ST using the format specified in FMT. Possible values of FMT
%   are 'mat' or 'nc'. We stronly recommend to use always file extensions 
%   .mat or .nc in STMFILE. Extensions other than .mat or .nc are not 
%   supported in data mode. 
%
%   STMWRITE(ST,STMFILE,'NODATA') writes the Space Time Matrix structure 
%   ST without obsData, auxData and stochData fields. This is useful
%   when used in combination with the data mode of this function.
%
%   STMWRITE(ST,STMFILE,'APPEND') is similar to 'NODATA', but should the 
%   STMFILE already exist, it will update all fields, except
%   the data fields. 
%
%   The 'NODATA' and 'APPEND' mode are intended to operate in combination 
%   with the data mode of the function.
%
%   In data mode the function takes as first argument a numeric Matlab
%   matrix. It assumes the STMFILE already exists (using the extension
%   .mat or .nc) and will only update one of the data fields.
%
%   STMWRITE(MAT,STMFILE,DATATYPE) writes the Matlab matrix MAT (Space Time
%   Matrix) with observation, auxiliary or stochastic data to the mat-
%   or NetCDF file STMFILE. DATATYPE is either 'OBSDATA', 'AUXDATA' or 
%   'STOCHDATA'. 
%
%   STMWRITE(MAT,STMFILE,DATATYPE,PNTIDX,EPOIDX,COMPIDX) writes a subset 
%   of DATATYPE to STMFILE. The Matlab matrix MAT (space Time
%   Matrix) contains the subset. The range of the indices is not
%   checked in this function. If indices are used, PNTIND and EPOIDX both
%   must be specified explicitly. If, COMPIDX is missing, a value of 1 (first
%   component) is assumed. This option is not supported for mat files 
%   below version V7.3.
%
%   The function can be used also with mat file objects STMOBJECT instead
%   of the matlab structure ST. It is fully compatible and uses the same
%   syntax. 
%
%   STMFILE=STMWRITE(STMOBJECT) returns only the file name STMFILE for the
%   STMOBJECT but does not do anything else. This is because the mat-file 
%   STMFILE is created at the time the matfile object is created using the 
%   STM function. STMFILE can be used in data mode to update the objects, 
%   but it is more efficient to pass the object itself, as shown below.
%
%   STMWRITE(MAT,STMOBJECT,DATATYPE) writes Matlab matrix MAT to the matfile 
%   object STMOBJECT as DATATYPE.
%
%   STMWRITE(MAT,STMOBJECT,DATATYPE,PNTIDX,EPOIDX,COMPIDX) writes a subset
%   of DATATYPE to the matfile object STMOBJECT.
%
%   Other calling syntaxes makes no sense for matfile objects. In these
%   cases a warning message is generated.
%
%   For creating Space Time Matrix datasets and more details on the dataset
%   stucture see the help of the STM function.
%
%   Examples:
%
%     st=stm(projectId,datasetId,techniqueId)
%     st.pntName={'a';'b';'c';'d'};        % fill the structure
%     st.epochDyear=2000:2019;
%     stmwrite(st,'test.mat','nodata')     % write the file (w/o data)
%     stmwrite(zeros(10,10),'test.mat','obsdata')
%     stmwrite(zeros(10,10),'test.mat','obsdata',16:25,3:12)
%     stmwrite(st,'test.mat','append')     % optionally update the structure
%
%   See also STM and STMREAD.
%
%   (c) Delft University of Technology, 2020.

% Created:  20 August 2020 by Hans van der Marel
% Modified: 25 August 2020 by Hans van der Marel
%            - Added data mode to the function
%            - Option to write subsets of the data
%            - Improved help text

% Check first input arguments and determine the main operation mode

if nargin <1, error('Need to have at least one input argument.'); end
if isobject(st)
   % Do nothing if first argument is an object, except returning the filename
   if isobject(st) && nargin > 1
      warning('Filename cannot be changed for STM objects, call function again with a single input argument, or in data mode.')
   end
   if nargout > 0
      varargout{1}=stmfile;
   end
   return
elseif not( isstruct(st) || isnumeric(st) )
   error('First argument must be a matlab structure, object or matrix.')
end

% Check the other arguments (depending of the mode)

if isstruct(st)

   % structure mode
   
   if nargin < 2
      datasetAttrib=st.datasetAttrib;
      stmfile=datasetAttrib.createdAs;
   else
      datasetAttrib=st.datasetAttrib;
      datasetAttrib.createdAs=stmfile;
      st.datasetAttrib=datasetAttrib;
   end
   if ~ischar(stmfile) && ~isstring(stmfile) 
      error('Second input argument (stmfile) must be a character string.')
   end
   
   if nargin < 3
      fmt=[];
      mode='DATA';
   else
      if any(strcmpi(mode,{'mat' 'nc'}))
         fmt=lower(mode);
      elseif any(strcmpi(mode,{'NODATA','APPEND'}))
         fmt=[];
      else
         error('Third argument must be a valid format, NODATA or APPEND')
      end
   end
   if isempty(fmt)      
      [~,~,ext]=fileparts(stmfile);
      fmt=lower(ext(2:end));
   end

   if nargin > 3
      error('Too many input arguments.')
   end
      
elseif isnumeric(st)
    
   % Data mode

   if nargin ~= 3 && nargin ~= 5 && nargin ~= 6
      error('Function requires 3, 5 or 6 arguments in data mode.')
   end
   
   if ischar(stmfile) || isstring(stmfile)
      [~,~,ext]=fileparts(stmfile);
      fmt=lower(ext(2:end));
      if ~any(strcmp(fmt,{'mat' 'nc'}))
         error('unsupported file format')
      end
      if ~exist(stmfile,'file')
         error(['Space time matrix file ' stmfile ' must already exist in datamode.'])
      end
   elseif isobject(stmfile)
      fmt='mat';
   else
      error('Second input argument (stmfile) must be a character string or object.')
   end  

   if ~any(strcmpi(mode,{'obsdata' 'auxdata','stochdata','nodata'}))
      error('unsupported data type')
   end
       
   dosubset=false;
   if nargin > 3 && nargin <= 6 && any(strcmpi({'OBSDATA' 'AUXDATA' 'STOCHDATA'},mode))
      if nargin < 5 
         error('PNTIDX and EPOIDX  must both be supplied')
      elseif nargin == 6
         compidx=varargin{3};
      else
         compidx=1;
      end
      pntidx=varargin{1};
      epochidx=varargin{2};
      dosubset=true;
   elseif nargin > 3 && nargin <= 6
      warning('PNTIDX, EPOIDX and COMPIDX only needed when MODE is a DATATYPE.') 
   elseif nargin > 6
      error('Too many input arguments.')  
   end

end


% Write the space time matrix to disk (for objects we are done...)

switch lower(fmt)
    case 'mat'
        % Read space time matrix from Matlab mat file
        switch upper(mode)
            case { 'DATA' }
               % Write space time matrix as Matlab mat file
               if exist(stmfile,'file')
                  disp(['Space time matrix dataset ' stmfile ' already exists, will be overwritten.'])
               end
               save(stmfile,'-v7.3','-struct','st')
            case { 'NODATA' }
               % Only save variables that do not match "Data" at end of name
               if exist(stmfile,'file')
                  disp(['Space time matrix dataset ' stmfile ' already exists, will be overwritten.'])
               end
               save(stmfile,'-v7.3','-struct','st','-regexp','...(?<!Data)$')
            case { 'APPEND' }
               % Only save variables that do not match "Data" at end of name
               if ~exist(stmfile,'file')
                  error(['Space time matrix dataset ' stmfile ' must already exist.'])
               end
               save(stmfile,'-append','-struct','st','-regexp','...(?<!Data)$')
            case { 'OBSDATA' 'AUXDATA' 'STOCHDATA' }
               % In these modes st is a numeric matrix, the matfile must already exist
               var=strrep(lower(mode),'data','Data');
               if dosubset
                  if isobject(stmfile)
                      stmfile.(var)(pntidx,epochidx,compidx)=st;
                  else
                      % create temporary object for stmfile with write properties
                      stmobj = matfile(stmfile,'Writable',true);
                      if numel(size(stmobj.(var))) == 2 && compidx == 1
                         stmobj.(var)(pntidx,epochidx)=st;
                      else
                         stmobj.(var)(pntidx,epochidx,compidx)=st;
                      end
                  end
               else
                  if isobject(stmfile)
                      stmfile.(var)=st;
                  else
                      % create temporary object for stmfile with write properties
                      stmobj = matfile(stmfile,'Writable',true);
                      stmobj.(var)=st;
                      % for matfiles older than v7.3 we can use save with -append
                  end
               end
            otherwise
               error(['Wrong mode ' mode ' (must be DATA, NODATA, OBSDATA, AUXDATA, STOCHDATA or OBJECT'])
        end
    case 'nc'
        % Read space time matrix from NetCDF file
        % stm=stmreadnetcdf(stmfile,mode)  % Internal function 
        error('NetCDF file format not yet supported')
    otherwise
        error('Unknown format')
end

% Output argument

if nargout > 0
    varargout{1}=stmfile;
end

end