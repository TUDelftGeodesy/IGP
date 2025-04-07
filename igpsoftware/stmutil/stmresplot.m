function stmresplot(outputfilename,varargin)
%stmresplot   Plot residuals and test statistics from stmintegrate.
%   STMRESPLOT(outputfilename) prints the test statistics, plots the point
%   test statistics, and plots the residuals, from stmintegrate.
%   The residuals are saved in the file OUTPUTFILENAME_res.mat.
%
%   STMRESPLOT(outputfilename,'option','value',...) with options
%
%      'doplots'      0 | 1                          Default 1
%      'visible'      'on'|'off'                     Default 'on'  *)
%      'saveplot'     'off'|'pdf'|'fig'|'pdf&fig'    Default 'off' *)
%      'background'   Map background                 Default 'igp-background.mat'
%      'unstablearea' Polygon defining unstable area Default 'igp-unstablearea.mat'
%
%  *) Not yet supported
%
%   See also STMINTEGRATE and STMRESIDUALS.
%
%  (c) Hans van der Marel, Delft University of Technology, 2021.

% Created:  02 March 2021 by Hans van der Marel
% Modified: 02 June 2022 by Hans van der Marel
%           - use pntIds as pntName if pntName variable does not exist
%           23 Oct 2023 by Hans van der Marel
%           - merged stmresplot2 with stmresplot (supporting options)
%           14 Nov 2023 by Hans van der Marel
%           - fixed bug (overide options after reading from residual file)

% Check arguments and process options

if nargin < 1
   error('Function must have at least one input argument')
end

% prepare outputId

[filepath,outputId] = fileparts(outputfilename);
if strcmpi(outputId(end-3:end),'_res')
   outputId=strrep(outputId,'_res','');
end

% read residuals from mat file

load(fullfile(filepath,[ outputId '_res.mat']));
if ~exist('pntName','var')
   pntName=pntIds;
end

% overide options loaded from the residuals file (this must be after
% reading the residuals file, do not move in front).

opt.doplots=1; 
opt.visible='on';
opt.saveplot='off';

if ~isfield(opt,'background')
   opt.background='igp-background.mat';     % map background
end
if ~isfield(opt,'unstablearea')
   opt.unstablearea='igp-unstablearea.mat'; % polygon defining unstable area
end

for k=1:2:numel(varargin)
   if isfield(opt,varargin{k})
      opt.(varargin{k})=varargin{k+1};
   else
      error(['Illegal option ' varargin{k} ])
   end
end

% call stmresiduals to do the plotting

if exist('pntCrd','var')
   stmresiduals(pntName, pntCrd, epochIds, dslocal, lsqstat1, opt, outputId)
else
   opt.doplots_mapcrd=true;
   stmresiduals(pntName, pntNeu, epochIds, dslocal, lsqstat1, opt, outputId)
end

end