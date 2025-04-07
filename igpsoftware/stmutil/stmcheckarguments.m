function [inputfilenames,outputfilename,opt]=stmcheckarguments(inputfilenames,outputfilename,defaultopt,varargin)
%stmcheckarguments   Check function arguments for integrated geodetic processing
%   [INPUTFILENAMES,OUTPUTFILENAME,OPT]=STMCHECKARGUMENTS(INPUTFILENAMES,
%   OUTPUTFILENAME,DEFAULTOPT,OPTIONS,FLAG) checks the INPUTFILENAMES and 
%   OUTPUTFILENAME, and sets the processing options. 
%
%   INPUTFILENAMES must be the name of an file containing the input filenames, 
%   or a cell array with the input filenames. On output INPUTFILENAMES is a 
%   cell array. If some of the input files do not exist the function 
%   returns an error. Files containing input filenames can have comment 
%   lines (lines starting with # or %).
%
%   OUTPUTFILENAME is the name of the output dataset. If the OUTPUTFILENAME
%   is empty on output, the OUTFILENAME exists and is up to date.
%
%   DEFAULTOPT contains the default processing options (set in the main 
%   program). OPTIONS is a structure, or cell array, with the options
%   to be changed. On output OPT is the combination of DEFAULTOPT and OPTIONS. 
%
%   FLAG is a keyword that affect the processing (opt.outputTarget) 
%   in case the target OUTPUTFILENAME is an existing file
%
%   flag       target OUTPUTFILENAME already exist
%   ---------  -----------------------------------------------------------
%   create     don't overwrite, stop processing 
%   update     overwrite if any of the input files is more recent, else stop processing
%   overwrite  overwrite target
%
%   If the target OUTPUTFILENAME does not yet exist this flag has no
%   effect. If no processing has to be done, OUTPUTFILENAME will be
%   empty.
%
%   This function should be called at the start of each module, directly
%   after defining default options.
%
%   Example:
%
%     function stmexample(inputfilenames,outputfile,varargin)
%     % Set Default options
%     defaultopt.a=1;
%     defaultopt.b=2;
%     % Check the input arguments
%     [inputfilenames,outputfile,opt]=stmcheckarguments(inputfilenames, ...
%           outputfile,defaultopt,varargin{:});
%     if isempty(outputfile)
%        return;
%     end
%     % Continue processing 
%     ...
%     end
%
%  The function stmexample is called as follows
%
%     stmexample(inputfilenames,outputfile,{'a',.5},'update')
%
%  with a new value for option 'a'. Actual processing will only be done
%  if one of the inputfiles is newer than the output file, or when
%  the options have been changed.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:  23 August 2020 by Hans van der Marel
% Modified: 26 August 2020 by Freek van Leijen, added cell to char
%               conversion
%           25 Sept 2020   by Freek van Leijen, small change to avoid
%               . indexing on functions (to enable use of older Matlab
%               versions) 
%            2 November 2020 by Hans van der Marel
%              - Added wildcards to inputfilenames
%           22 October 2021 by Hans van der Marel
%              - Ability to comment out lines in files containing
%                file names, lines starting with # or % are treated 
%                as comments.
%              - Improved print output of structures
%           02 November 2021 by Freek van Leijen
%              - corrected a nested loop with the same loop-counter
%              - added option output for a nested cell array
%            4 November 2021 by Hans van der Marel
%              - repaired printing of options
%           26 October 2022 by Hans van der Marel
%              - output for structure options by calling stmprintopt

% Check the input arguments

if nargin < 3
    error('This function expects at least three input arguments.')
end

% Check the options and if necessary overwrite the default values

opt=defaultopt;
updatemode='create';
for l=1:numel(varargin)
   if isstruct(varargin{l}) || iscell(varargin{l})
      if iscell(varargin{l})
         options=cell2struct(varargin{l}(2:2:end),varargin{l}(1:2:end),2);
      else
         options=varargin{l};
      end
      modifiedoptions=fieldnames(options);
      for k=1:numel(modifiedoptions)
         modifyoption=modifiedoptions{k};
         if isfield(opt,modifyoption)
            fprintf('Option %s changed: ',modifyoption)
            % Print default value (opt)
            stmprintopt(opt.(modifyoption))
            fprintf(' -> ')
%             if isempty(opt.(modifyoption))
%                fprintf(' <empty> -> ')
%             elseif iscell(opt.(modifyoption))
%                tmp=opt.(modifyoption);
%                for m=1:numel(opt.(modifyoption))
%                    if iscell(opt.(modifyoption){m})
%                       tmp{m}=['\n  - ' strjoin(opt.(modifyoption){m})];
%                    end
%                end
%                fprintf('%s -> ',strjoin(tmp)) 
%             elseif isstruct(opt.(modifyoption))
%                fprintf('\n')
%                %disp(opt.(modifyoption)) 
%                if numel(opt.(modifyoption)) < 1
%                    fprintf(' <empty>')
%                end
%                for m=1:numel(opt.(modifyoption))
%                   disp(opt.(modifyoption)(m))
%                end
%                fprintf(' -> ')
%             elseif ischar(opt.(modifyoption)) || isstring(opt.(modifyoption))
%                fprintf('%s -> ', opt.(modifyoption))
%             else
%                fprintf('%s -> ', mat2str(opt.(modifyoption)))
%             end
            % Print modified value (options)  
            if isstruct(options.(modifyoption))
               stmprintopt(options.(modifyoption),'   ')
            else
               stmprintopt(options.(modifyoption))
            end
            fprintf('\n')
%             if isempty(options.(modifyoption))
%                fprintf(' <empty>\n')
%             elseif iscell(options.(modifyoption))
%                tmp=options.(modifyoption);
%                for m=1:numel(options.(modifyoption))
%                    if iscell(options.(modifyoption){m})
%                       tmp{m}=['\n  - ' strjoin(options.(modifyoption){m})];
%                    end
%                end
%                fprintf('%s\n',strjoin(tmp)) 
%             elseif isstruct(options.(modifyoption))
%                fprintf('\n')
%                for m=1:numel(options.(modifyoption))
%                   disp(options.(modifyoption)(m))
%                end
%             elseif ischar(options.(modifyoption)) || isstring(options.(modifyoption))
%                fprintf('%s\n', options.(modifyoption))
%             else
%                fprintf('%s\n', mat2str(options.(modifyoption)))
%             end                
            opt.(modifyoption)=options.(modifyoption);
         else
            warning(['Invalid option ' modifyoption ', ignore...'])
         end
      end
   elseif ~isempty(varargin{l})
      switch lower(varargin{l})
         case {'create' 'update' 'overwrite' }
            updatemode=lower(varargin{l});
         otherwise
            error(['Invalid commandline option ' varargin{l} ])
      end
   end
end
opt.outputTarget=lower(updatemode);

% Check if opt.inputDir exist, if not, make it empty
if ~isfield(opt,'inputDir')
   opt.inputDir='';
end

% Check the input dataset names

if ischar(inputfilenames) || isstring(inputfilenames)
   d=dir(fullfile(opt.inputDir,inputfilenames));
   if numel(d) > 1
      inputfilenames={d.name};
      opt.inputDir=char(unique({d.folder}));
   elseif numel(d) == 1
      filewithinputfiles=inputfilenames;
      fid=fopen(filewithinputfiles);
      k=0;
      inputfilenames=[];
      while ~feof(fid)
         line=strtrim(fgetl(fid));
         if isempty(line) || line(1) == '#' || line(1) == '%' || line(1) == '[' 
             continue;
         end
         k=k+1;
         inputfilenames{k}=line;
      end
      fclose(fid);
   else
      error(['The file ' inputfilenames ' with the inputfile names does not exist.'])
   end
%    if ~exist(filewithinputfiles,'file')
%       error(['The file ' filewithinputfiles ' with the inputfile names does not exist.'])
%    end  
%    fid=fopen(filewithinputfiles);
%    k=0;
%    inputfilenames=[];
%    while ~feof(fid)
%       k=k+1;
%       inputfilenames{k}=fgetl(fid);
%    end
%    fclose(fid);
elseif ~iscell(inputfilenames)
   error('Inputfilenames must be a cell array with inputfile names or a string with the name of a file containing the inputfile names.')
end

% Check that the input datasets exist 

haveerror=false;
for k=1:numel(inputfilenames)
   if ~exist(fullfile(opt.inputDir,char(inputfilenames{k})),'file')
      disp(['The input dataset ' char(inputfilenames{k}) ' does not exist.'])
      haveerror=true;
   end  
end
if haveerror
    error('One or more of the input datasets do not exist.')
end

% Check if the output dataset exists and may be overwritten
%
% outputTarget   target does not exist   target already exist
% -------------  ---------------------   ----------------------
% create         create                  don't overwrite, stop processing 
% update         create                  overwrite if any of the input files is more recent, if not stop processing
% overwrite      create                  overwrite

if exist(outputfilename,'file')
   switch lower(updatemode) 
       case 'update'
          outputDir=dir(outputfilename);
          moddateoutput=outputDir.datenum;
          havenewer=0;
          for k=1:numel(inputfilenames)
             inputDir=dir(inputfilenames{k});
             if inputDir.datenum > moddateoutput
                havenewer=havenewer+1;
             end
          end
          % also include a check if options have changed, by comparing with
          % options stored in the output dataset (TO BE IMPLEMENTED)
          if havenewer > 0
             disp(['There are ' num2str(havenewer) ' input datasets more recent that the output dataset, the output dataset ' outputfilename ' will be updated.'])
          else
             disp(['The output dataset ' outputfilename ' is up to date, no processing needed, use outputTarget=overwrite to overwrite the output dataset.'])
             % make outputfilename empty to signal to stop the processing
             outputfilename=[];
             return
          end
       case 'overwrite'
          disp(['The output dataset ' outputfilename ' exists and will be overwritten.'])
       case 'create'
          error(['The output dataset ' outputfilename ' exists, rename, or set the outputTarget setting, and restart the processing.'])
       otherwise
          error(['Incorrect setting ' updatemode ' for outputTarget (must be create, update or overwrite), abort processing.'])
   end  
end

end
