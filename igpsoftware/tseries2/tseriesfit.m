function [tseries,hf]=tseriesfit(station,varargin)
%TSERIESFIT  Fit time series components.
%   TSERIES=TSERIESFIT(TSERIESIN,OPT) fits the models specified in
%   the structure OPTIN to the time series data TSERIESIN. Output is a
%   new time series structure TSERIES with the fitted components.
%
%   TSERIES=TSERIESFIT(STATION,OPT) fits the time series data for station 
%   STATION, with the time series data stored in the Matlab mat file
%   <STATION>.mat.
%
%   TSERIES=TSERIESFIT(STATION|TSERIESIN,OPTNAME,OPTVALUE,...) uses option
%   pairs instead of structure.
%
%   TSERIES=TSERIESFIT(...,ID) also saves TSERIES to a file
%   <STATION><ID>.mat.
%
%   The models are specified in the Matlab structure OPT:
%
%     opt.trendmodel=1          %  trend model (default 'auto'), inputs
%                               %    #degree       polynominal degree
%                               %    'spline'      spline fit
%                               %    'auto'        linear < 2years or spline fit
%     opt.harmonic=[]           %  harmonic periods in decimal years (default [1 1/2.08]) 
%                               %    [1 1/2.08]    annual and semi-annual
%                               %    [1 347/365 0.5 14.2/365 ]
%     opt.meteo=[]              %  meteo data structure or meteo data filename
%     opt.rstack=[]             %  residual stack structure or residual stack file name
%     opt.minstack=2;           %  minimum number of contributing series to stack
%     opt.cm=[]                 %  common mode structure or common mode filename
%     opt.breaks=[]             %  insert breakpoints in decimal years (default [])
%
%     opt.maxsigma=[Inf Inf Inf]%  maximum a-priori sigma in [mm] (default Inf)
%                               %    [ 2 2 7 ]    remove everything with st.dev. > 2 mm horz and 7 mm vertical
%                               %    [ 1.5 1.5 5 idem with 1.5 and 5 mm
%     opt.maxresid=[Inf Inf Inf]%  maximum residual in [mm] for outlier detection
%                               %    [ 10 10 20 ]
%     opt.maxiter=2             %  maximum number of iterations for outlier removal
%
%     opt.doplots=0             %  plotting (0:none, 1: basic, 2:extensive, if negative, visible is off
%     opt.saveplot=false        %  save the plots to disk
%     opt.txtfile=false         %  save txt summary file to disk
%
%   See also tseriesanalysis.
%
%   (c) Hans van der Marel, Delft University of Technology, 2014-2019.

% Default options

opt.trendmodel='auto';       % auto detect linear or spline fit by default
opt.harmonic=[1 1/2.08 ];    % estimate annual and semi-annual by default
opt.meteo=[];                % meteo data structure or meteo data filename
opt.rstack=[];               % residual stack structure or residual stack file name
opt.minstack=2;              % minimum number of contributing series to stack
opt.cm=[];                   % common mode structure or common mode filename
opt.breaks=[];               % do not insert break points

opt.maxsigma=[Inf Inf Inf];  % do not remove points with bad a-priori st.dev.  
opt.maxresid=[Inf Inf Inf];  % maximum residual in [mm] for outlier detection
opt.maxiter=2;               % maximum number of iterations for outlier removal

opt.debias=false;

opt.doplots=0;               % by default we don't do plotting, 1 is basic plotting, 2 extensive
opt.saveplots=false;         % by default plots are not saved 
opt.txtfile=false;           % by default, no txt files are saved

% Check the input arguments

if ~ischar(station) && ~isstruct(station)
   error('First input argument must be a tseries structure or station name.')
end
if ( nargin < 2 ) || (nargin <= 3 && isstruct(varargin{1})) 
   if nargin < 2
     optin=struct();
   else
     optin=varargin{1};
   end
   if nargin < 3
     savename=[];
   else
     savename=varargin{2};
   end
else
   for k=1:2:numel(varargin)-1
      optin.(varargin{k})=varargin{k+1};
   end
   if numel(varargin) > 2*floor(numel(varargin)/2)
     savename=varargin{end};
   else
     savename=[];
   end
end
if isempty(savename)
  savedata=false;
  savename='_fit';
else
  savedata=true;  
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

% Get the data

if ischar(station)
   load([ station '.mat' ]);
elseif isstruct(station)
   tseries=station;
   station=tseries.station;
end

epoch=tseries.epoch;
if isfield(tseries,'year')
  year=tseries.year;
else
  year=date2dyear(epoch);
end
neu=tseries.neu;
sneu=tseries.sneu;

if isfield(tseries,'id')
  id=tseries.id;
else
  id='';
end
 
% Get antenna and receiver names and other events

%  EVENTS is a structure array with at least the following fields
%    events(k).year   start of the k'th event in decimal years
%    events(k).type   type of event: REC=receiver change, ANT=antenna
%                     change, or others.
%    events(k).name   name of the receiver or antenna (only for REC or ANT)

if isfield(tseries,'events')
   events=tseries.events;
else
   events=[];
end
if iscellstr(opt.breaks)
   if ~isempty(events)
     eventype={ events(:).type }';
     idx=ismember(eventype,opt.breaks);
     breaks=[ events(idx).year ]';
     brkname=[ events(idx).name ]';
     itmp=find(breaks < year(1)+2/365 | breaks > year(end)-2/365);
     breaks(itmp)=[];    
     brkname(itmp,:)=[];
     [breaks,itmp]=sort(breaks);
     brkname=brkname(itmp,:);
   else
     breaks=[];
   end
else
   breaks=opt.breaks;
   % add break to event list
   kk=numel(events);
   for k=1:numel(breaks)
      events(kk+k).year=breaks(k);
      events(kk+k).type='STEP';
      events(kk+k).name='';
   end      
end

% Get nominal position and nominal velocity

plh0=tseries.plh0;          % latitude and longitude are in degrees
if isfield(tseries,'vneu0')
  vneu0 =tseries.vneu0;
else
  vneu0=[0 0 0];
end
plh=[ plh0(1:2)*pi/180 plh0(3)];

% Filter data points on a-priori sigma

maxsigma=opt.maxsigma./1000;

idx=find(any([ sneu(:,1) > maxsigma(1) sneu(:,2) > maxsigma(2) sneu(:,3) > maxsigma(3) ],2));

if numel(idx) > 0

  fprintf('Remove data points with large St.Dev. (sN > %.1f , sE > %.1f  sU %.1f [mm]):\n\n',maxsigma*1000)
  fprintf(' Epoch              North [mm]  East [mm]    Up [mm]    sN [mm]  sE [mm]  sU [mm]\n')
  for k=1:length(idx)
    fprintf('%6d (%8.3f)   %10.1f %10.1f %10.1f  %8.1f %8.1f %8.1f\n',idx(k),year(idx(k)),neu(idx(k),:)*1000,sneu(idx(k),:)*1000);
  end

  epoch(idx)=[];
  year(idx)=[];
  neu(idx,:)=[];
  sneu(idx,:)=[];

end

% Correct input timeseries with nominal velocity

% dneu=neu;                 % dneu is the series with corrections

dneu=neu-(year-mean(year))*vneu0;

% Optionally correct input timeseries with residual stack

if ~isempty(opt.rstack)
  disp('correct NEU with residual stack...')
  if isstruct(opt.rstack)
     rstack=opt.rstack;
  elseif ischar(opt.rstack)
     load(opt.rstack);
  else
     error('Residual stack option must be structure or a filename.')
  end
  ghour=round((epoch-datenum([2013 1 1]))*24);
  [~,idx]=ismember(ghour,rstack.rhour);
  if any(idx==0)
    %ghour(idx==0)
    disp([ num2str(numel(find(idx==0))) ' data points are not in the residual stack, remove...'])
  end
  idx2=(idx~=0);
  idx3=find(any(rstack.rcount(idx(idx2)) < opt.minstack,2));
  if numel(idx3) > 0
    disp([ num2str(numel(idx3)) ' data points have a residual stack based on less than ' num2str(opt.minstack) ' series  , remove...'])
    idx2(idx3)=false;
  end
  epoch=epoch(idx2);
  year=year(idx2);
  neu=neu(idx2,:);
  dneu=dneu(idx2,:)-rstack.rmean(idx(idx2),:);
  sneu=sneu(idx2,:);
end

% Load meteo data

metdata=[];  
if ~isempty(opt.meteo)
  disp('Get meteo data...')
  if isstruct(opt.meteo)
     metstruct=opt.meteo;
  elseif ischar(opt.meteo)
     metstruct=load(opt.meteo);
  else
     error('Meteo data option must be structure or a filename.')
  end
  d2=datestr(epoch,'yyyymmdd');
  [~,k2]=ismember(d2,metstruct.YYYYMMDD,'rows');
  metdata={ metstruct.Pday(k2) metstruct.Tday(k2) };
end

% Correct for common mode

if ~isempty(opt.cm)
  disp('correct NEU with common mode...')
  if isstruct(opt.cm)
     cm=opt.cm;
  elseif ischar(opt.cm)
     load(opt.cm);
  else
     error('Common mode option must be structure or a filename.')
  end
  if isempty(opt.meteo)
     error('In order to use common mode corrections you must load meteo data.')
  end
  ycm=tseriescmeval(cm,year,metdata);
  dneu=dneu-ycm;
end


% Analyze the time series with outlier detection (data cleaning)

for l=1:opt.maxiter

  % Apply corrections for timeseries steps

  if opt.debias

     dneu2=dneu;
     [dneu2(:,1),ibrk,jump1]=debias(year,dneu2(:,1),breaks);
     dneu2(:,1)=dneu2(:,1)-mean(dneu2(:,1));
     [dneu2(:,2),ibrk,jump2]=debias(year,dneu2(:,2),breaks);
     dneu2(:,2)=dneu2(:,2)-mean(dneu2(:,2));
     [dneu2(:,3),ibrk,jump3]=debias(year,dneu2(:,3),breaks);
     dneu2(:,3)=dneu2(:,3)-mean(dneu2(:,3));
     breaks=year(ibrk(2:end-1));
       
     fprintf('\nTime series shifts (debias):\n\n')
     fprintf(' Epoch (year)       North [mm]  East [mm]    Up [mm]\n')
     for k=1:length(ibrk)-1
       fprintf('%6d (%8.3f)   %10.1f %10.1f %10.1f   %s\n',ibrk(k),year(ibrk(k)),jump1(k)*1000,jump2(k)*1000,jump3(k)*1000,brkname{k});
     end
     fprintf('%6d (%8.3f)\n',ibrk(end)-1,year(ibrk(end)-1))

     dneu=dneu2;
  end
  
  % Detrend and estimate individual components

  tseries=tseriesanalysis(year,dneu,sneu,opt.trendmodel,metdata,opt.harmonic,breaks);

  tseries.station=station;
  tseries.id=id;
  tseries.plh0=plh0;
  tseries.vneu0=vneu0;
  tseries.events=events;

  neuresidual=tseries.neuresidual;
  neunojumps=tseries.neu-tseries.neujumps;
  neufit=tseries.neufit-tseries.neujumps;
  neutrend=tseries.neutrend;
  neutrendldn=tseries.neutrend+tseries.neuatmld;

  % Remove bad data points

  if ischar(opt.maxresid) 
    idx=tseriesoutlier(tseries,'item','residual','crit',5,'stepcrit',5,'verbose',1);
    if (isempty(idx)) 
      break;
    end
  else
      
    maxdelta=opt.maxresid./1000;

    idx=find(any([ abs(neuresidual(:,1)) > maxdelta(1) abs(neuresidual(:,2)) > maxdelta(2) abs(neuresidual(:,3)) > maxdelta(3) ],2));

    if (isempty(idx)) 
      fprintf('\nNo outliers found (dN > %.1f , dE > %.1f  dU > %.1f [mm]), continue...\n\n',maxdelta*1000)
      break;
    end

    fprintf('\nRemove outliers (dN > %.1f , dE > %.1f  dU > %.1f [mm]):\n\n',maxdelta*1000)
    fprintf(' Epoch              North [mm]  East [mm]    Up [mm]    sN [mm]  sE [mm]  sU [mm]\n')
    for k=1:length(idx)
      fprintf('%6d (%8.3f)   %10.1f %10.1f %10.1f  %8.1f %8.1f %8.1f\n',idx(k),year(idx(k)),dneu(idx(k),:)*1000,sneu(idx(k),:)*1000);
    end
  end
    
  year(idx)=[];
  neu(idx,:)=[];
  dneu(idx,:)=[];
  sneu(idx,:)=[];
  neuresidual(idx,:)=[];
  neunojumps(idx,:)=[];
  neufit(idx,:)=[];
  neutrend(idx,:)=[];
  neutrendldn(idx,:)=[];
  if ~isempty(metdata)
    P=metdata{1};
    T=metdata{2};
    T(idx)=[];
    P(idx)=[];
    metdata={ P T };
  end

end

% Compute scaling factors for standard deviations

empstd=std(diff(dneu))./sqrt(2)*1000;
meanstd=mean(sneu)*1000;

f=empstd./meanstd;

fprintf('\n                    North [mm]  East [mm]    Up [mm]\n')
fprintf('Emperical St.Dev.:  %10.3f %10.3f %10.3f\n',empstd);
fprintf('Formal St.Dev.:     %10.3f %10.3f %10.3f\n',meanstd);
fprintf('Factor (Estimated): %10.3f %10.3f %10.3f\n',f);

%f=[ 1 1 1];
%fprintf('Factor (Applied):   %10.3f %10.3f %10.3f\n\n',f);

if opt.txtfile
  fid=fopen([ station savename '.txt' ],'w');
else
  fid=1;
end
fprintf(fid,'\n\n           Vframe  Vsite   +/- [mm/y]');
fmt='';
for k=1:length(tseries.harmonic)
   fprintf(fid,'%3dd   ',round(tseries.harmonic(k)*365));
   fmt=[ fmt ' %6.2f'];
end
fprintf(fid,' StdR   StdE   StdF [mm] \n');
fprintf(fid,[ '%s  Lat  %6.2f %6.2f %6.2f   ' fmt '  %6.2f %6.2f %6.2f     %7.2f-%7.2f\n%s  Lon  %6.2f %6.2f %6.2f   ' fmt '  %6.2f %6.2f %6.2f     %7.2f-%7.2f\n%s  Hgt  %6.2f %6.2f %6.2f   ' fmt '  %6.2f %6.2f %6.2f     %7.2f-%7.2f\n\n'],...
    station,vneu0(1)*1000,tseries.vel(1)*1000,tseries.svel(1)*1000,tseries.amplitude(1,:)*1000,tseries.rms(1)*1000,empstd(1),meanstd(1),tseries.tfirst,tseries.tlast, ...
    station,vneu0(2)*1000,tseries.vel(2)*1000,tseries.svel(2)*1000,tseries.amplitude(2,:)*1000,tseries.rms(2)*1000,empstd(2),meanstd(2),tseries.tfirst,tseries.tlast, ... 
    station,vneu0(3)*1000,tseries.vel(3)*1000,tseries.svel(3)*1000,tseries.amplitude(3,:)*1000,tseries.rms(3)*1000,empstd(3),meanstd(3),tseries.tfirst,tseries.tlast);
fprintf(fid,'%s  PLH  %12.8f %12.8f %9.4f    (%s )\n',station,plh(1:2)*180/pi,plh(3),plh2str(plh));
if opt.txtfile
  fclose(fid);
  type([ station savename '.txt']);
end

% Save time series results

if savedata
  disp(['save tseries fit to file ' station savename '.mat'])
  save([ station savename '.mat'],'tseries');
end

% Plot the results

hf=tseriesplot(tseries,'doplots',opt.doplots,'saveplots',opt.saveplots,'plotdir','plots','empstdev',false);

end


function [xout,ibrk,jump]=debias(year,xin,breaks)

for i=1:length(breaks)
   idx=find(year > breaks(i));
   ibrk(i)=idx(1);
end
ibrk(length(breaks)+1)=length(year)+1;
if ibrk(1) > 1
    ibrk=[1 ibrk];
end

jump(1)=0;
for i=2:length(ibrk)-1
   i1=[max([ibrk(i)-31 ibrk(i-1)]) ibrk(i)-1 ];
   i2=[ibrk(i) min([ibrk(i)+30 ibrk(i+1)-1 ])];
   jump(i)=mean(xin(i2(1):i2(2)))-mean(xin(i1(1):i1(2)));
   if isnan(jump(i)); jump(i)=0; end
end
tjump=0;
xout=xin;
for i=1:length(jump)
   i1=ibrk(i);
   i2=ibrk(i+1)-1;
   tjump=tjump+jump(i);
   xout(i1:i2)=xout(i1:i2)-tjump;
end

end

