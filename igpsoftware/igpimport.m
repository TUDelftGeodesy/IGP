function igpimport(toolbox,environment)
%IGPIMPORT Add IGP toolbox to the Matlab path. 
%   IGPIMPORT(TOOLBOX) adds the directory with TOOLBOX to the 
%   Matlab path. TOOLBOX is the toolbox name, the actual directory path
%   is defined in the file igptoolbox.cfg.
%
%   IGPIMPORT(TOOLBOX,ENVIROMENT) adds the directory with TOOLBOX to the 
%   Matlab path, but now using the directory paths defined in
%   ENVIRONMENT.cfg file.
%
%   IGPIMPORT() lists the available toolboxes.
%
%   (c) Hans van der Marel, Delft University of Technology, 2021.

%   Created:    15 Sep 2021 by Hans van der Marel
%   Modified:   

% Find out the path to this file

mfilename;
thisfile=which(mfilename);
rootdir=fileparts(thisfile);

% Name of the file with path definitions

if nargin < 2
    environment='igptoolbox';
end
envfile = fullfile(rootdir,[ environment '.cfg']);
if ~exist(envfile','file')
    error(['IGP toolbox configuration file ' envfile ' does not exist.' ])
end

toolboxdir = '';

% Read file with path definitions 

fid = fopen(envfile);
while ~feof(fid)
    line = fgetl(fid);
    if nargin < 1
        disp(line)
    else
        keyvalue = split(line,'=');
        if strcmpi(toolbox,strip(keyvalue{1}))
            toolboxdir=strip(keyvalue{2});
        end
    end    
end
fclose(fid);

if nargin < 1
    return
end

% Add toolboxdir to the Matlab path

if strcmpi(toolboxdir,'')
    error(['IGP tooolbox ' toolbox ' not in configuration file ' envfile ])
end
if ~exist(toolboxdir,'dir')
    error(['IGP toolbox directory ' toolboxdir ' does not exist.' ])
end

addpath(toolboxdir)
fprintf('Added toolbox %s (%s) to the Matlab path\n',toolbox,toolboxdir);

end
