function varargout=tseriesoutlier(ts,varargin)
%TSERIESOUTLIER  Robust outlier and step detection for time series objects.
%   TSERIESOUTLIER(TS) does outlier and step detection on the timeseries 
%   TS, using a robust moving average filter implemented in RMAFILT. TS 
%   is a time series structure or structure array.
%
%   [IDXOUTLIER,IDXSTEP,STEPS]=TSERIESOUTLIER(TS) outputs the index array
%   with detected outliers IDXOUTLIER and steps IDXSTEP, and an approximation
%   of the step size. The output is only produced when TS is a structure
%   or single element of an structure array.
%
%   [...]=TSERIESOUTLIER(...,OPTION, VALUE) passes OPTION/VALUE pairs, with 
%
%     item      cell array with items to be used as input for outlier 
%               detection. The default is 'residual' in case of a fitted
%               stucture, or else 'neu' in case of raw-data.
%     verbose   verbose level (default 1)
%                  0    no output
%                  1    print of detected outliers and steps
%                  2    more extensive output
%                  3    internal output from RMAFILT
%                  4    plotting by RMAFILT
%     window    window length for filter (default 21)
%     crit      outlier detection limit (default 5)
%     stepcrit  step detection limit (default 5), if 0 or NaN, no step 
%               detection is done.
%
%   See also RMAFILT.
%
%   (c) Hans van der Marel, Delft University of Technology, 2019.

%   Created:    28 January 2019 by Hans van der Marel
%   Modified:    5 September 2019 by Hans van der Marel
%                  - corrected the options

% Default options

opt.item='residual';
opt.verbose=1;
opt.window=21;
opt.crit=5;
opt.stepcrit=5;

% Check the input arguments

if ~isstruct(ts)
   error('First input argument must be a tseries structure.')
end
if isfield(ts,'neuresidual') 
  opt.item='residual';
else
  opt.item='neu';
end
if nargin < 2 || isstruct(varargin{1}) 
   if nargin < 2
     optin=struct();
   else
     optin=varargin{1};
   end
else
   for k=1:2:numel(varargin)-1
      optin.(varargin{k})=varargin{k+1};
   end
end

% Overide default options with the input options

options=fieldnames(optin);
if ~all(ismember(fieldnames(optin),fieldnames(opt)))
    disp('Warning: illegal field in input option:')
    options(~ismember(fieldnames(optin),fieldnames(opt)))
end
for k=1:numel(options)
   field=options{k};
   opt.(field)=optin.(field);
end

% select item for outlier detection

opt.item=cellstr(opt.item);
for k=1:length(opt.item)
  if ismember(opt.item,{'neu' 'raw' 'nostep'})
    opt.item{k}='neu';
  else
    opt.item{k}= [ 'neu' opt.item{k} ];
  end
end

remjump=false;
if ismember(opt.item,{'nostep' 'fit'})
  remjump=true;
end

% outlier detection

for k=1:numel(ts)

  if opt.verbose > 1
     fprintf('\n\nOutlier detection for %s\n\n',ts(k).station); 
     prtevents(ts(k));
  end

  data=ts(k).(opt.item{1});
  if remjump, data=data-ts.neujumps; end
  for l=2:length(opt.item)
    data=data+ts(k).(opt.item{l});
  end 

  [ym1,sy1,idxoutlier1,idxstep1,steps1,statoutlier1]=rmafilt(data(:,1)*1000,opt.window,'crit',opt.crit,'stepcrit',opt.stepcrit,'verbose',opt.verbose-2);
  [ym2,sy2,idxoutlier2,idxstep2,steps2,statoutlier2]=rmafilt(data(:,2)*1000,opt.window,'crit',opt.crit,'stepcrit',opt.stepcrit,'verbose',opt.verbose-2);
  [ym3,sy3,idxoutlier3,idxstep3,steps3,statoutlier3]=rmafilt(data(:,3)*1000,opt.window,'crit',opt.crit,'stepcrit',opt.stepcrit,'verbose',opt.verbose-2);

  idxoutlier=unique([idxoutlier1 ; idxoutlier2 ; idxoutlier3]);
  statoutlier=nan(size(idxoutlier,1),3,3);
  [~,ia,ib] = intersect(idxoutlier,idxoutlier1);
  statoutlier(ia,:,1)=statoutlier1(ib,:);
  [~,ia,ib] = intersect(idxoutlier,idxoutlier2);
  statoutlier(ia,:,2)=statoutlier2(ib,:);
  [~,ia,ib] = intersect(idxoutlier,idxoutlier3);
  statoutlier(ia,:,3)=statoutlier3(ib,:);
  
  idxstep=unique([idxstep1 ; idxstep2 ; idxstep3]);
  steps=nan(size(idxstep,1),3);
  if (~isempty(idxstep)) 
    steps=nan(size(idxstep,1),3);
    [~,ia,ib] = intersect(idxstep,idxstep1);
    steps(ia,1)=steps1(ib);
    [~,ia,ib] = intersect(idxstep,idxstep2);
    steps(ia,2)=steps2(ib);
    [~,ia,ib] = intersect(idxstep,idxstep3);
    steps(ia,3)=steps3(ib);
    if opt.verbose > 0
      fprintf('Steps (%d) found in %s\n\n',numel(idxstep),ts(k).station);
      fprintf('Epoch           North [mm] East [mm]   Up [mm]\n')
      for l=1:length(idxstep)
        fprintf('%5d(%8.3f) %10.1f %10.1f %10.1f\n',idxstep(l),ts(k).year(idxstep(l)),steps(l,:));
      end
    end
  end
  
  if (isempty(idxoutlier))
    if opt.verbose > 0
      fprintf('\nNo outliers found in %s\n\n',ts(k).station);
    end
    continue;
  end

  if opt.verbose > 0
    fprintf('\nOutliers (%d) found in %s\n\n',length(idxoutlier),ts(k).station)
    fprintf('Epoch           North [mm]________________ East [mm]_________________ Up [mm]___________________\n')
    for l=1:length(idxoutlier)
      line=sprintf('%5d(%8.3f) %5.1f%6.1f%6.1f%6.2f%3d %5.1f%6.1f%6.1f%6.2f%3d %5.1f%6.1f%6.1f%6.2f%3d\n',idxoutlier(l),ts(k).year(idxoutlier(l)),...
        ts(k).sneu(idxoutlier(l),1)*1000,data(idxoutlier(l),1)*1000,statoutlier(l,:,1), ...
        ts(k).sneu(idxoutlier(l),2)*1000,data(idxoutlier(l),2)*1000,statoutlier(l,:,2), ...
        ts(k).sneu(idxoutlier(l),3)*1000,data(idxoutlier(l),3)*1000,statoutlier(l,:,3));
      fprintf('%s',strrep(line,'NaN','   '));
    end
  end

end

if numel(ts) == 1
   varargout={ idxoutlier , idxstep ,steps };
end

end