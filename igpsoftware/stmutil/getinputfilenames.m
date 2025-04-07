function inputfilenames=getinputfilenames(filewithinputfiles,section)
%getinputfilenames   Read input files names from file.
%   INPUTFILENAMES=GETINPUTFILENAMES(FILEWITHINPUTFILES) reads the input 
%   file names from the file FILEWITHINPUTFILES and returns the file
%   names in the cell array INPUTFILENAMES. Lines starting with starting 
%   with a "#" or "%" character are treated as comment lines.
%
%   INPUTFILENAMES=GETINPUTFILENAMES(FILEWITHINPUTFILES,SECTION) only 
%   returns filenames inside a section starting with [SECTION].
%
%   Example file (example_file.txt)
%
%      # Example of comment line
%      [campaign]
%      ../../igpdata/nam_levelling/stm/gron_levelling_flaggedOutliers_v2.3.02_20191030.mat
%      [gnss]
%      ../../igpdata/nam_gnss/06gps_nam_202008.mat
%      [insar]
%      ../../igpdata/nam_insar/stm/U05_NAM_GTZH_U05_deformation.mat
%      ../../igpdata/nam_insar/stm/noord_nl_s1_asc_gaussian90_deformation.mat
%      ../../igpdata/nam_insar/stm/noord_nl_s1_dsc_gaussian90_deformation.mat
%
%    The sections are optional.
%
%    Examples
%
%      all_files=getinputfilenames('example_file.txt')
%      insar_files=getinputfilenames('example_file.txt','insar')
%
%  (c) Hans van der Marel, Delft University of Technology, 2021.

% Created:  22 October 2021 by Hans van der Marel
% Modified: 

if nargin < 2
   section=[];
end

inputfilenames=[];
currentsection=[];

k=0;
fid=fopen(filewithinputfiles);
while ~feof(fid)
  line=strtrim(fgetl(fid));
  if isempty(line) || line(1) == '#' || line(1) == '%'
     continue;
  end
  if line(1) == '[' && line(end) == ']'
      currentsection = line(2:end-1);
      continue;
  end
  if isempty(section) || strcmpi(currentsection,section) 
     k=k+1;
     inputfilenames{k}=line;
  end
end
fclose(fid);

inputfilenames=inputfilenames(:);

end