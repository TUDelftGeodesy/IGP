function [ym,sy,idxoutlier,idxstep,steps,statoutlier]=rmafilt(y,nwindow,varargin)
%RMAFILT  Robust moving average filter with outlier and step detection.
%   [YM,SY,IDXOUTLIER,IDXSTEP,STEPS]=RMAFILT(Y,NWINDOW) computes a robust 
%   moving average YM and standard deviation SY from the input data series
%   Y. NWINDOW is the window length for the moving average. IDXOUTLIER 
%   is an index vector with the outliers Y(IDXOUTLIER) and IDXSTEP is
%   an index vector with detected steps STEPS.
%
%   [...]=RMAFILT(...,OPTION, VALUE) passes OPTION/VALUE pairs, with 
%
%     verbose   verbose level (default 0)
%                  0 no output
%                  1 printed output
%                  2 with plots
%     crit      outlier detection limit (default 3.29)
%     stepcrit  step detection limit (default 1.5), if 0 or NaN, no step 
%               detection is done.
%
%   Examples:
%
%   (c) Hans van der Marel, Delft University of Technology, 2019.

%   Created:    28 January 2019 by Hans van der Marel
%   Modified:    5 September 2019 by Hans van der Marel
%                  - corrected the options

%% Dummy code section for function development (executed only in workspace from code editor)

if isempty(mfilename)

  % Set Length of the series, simulated steps and outliers 

  n=1000;
  sigma0=1.5;
            %  epoch magnitude
  simsteps    = [ 300     3 ; ...
                  550    -5 ; ...
                  610     3 ; ...
                  630    -5 ];
  simoutliers = [ 200     5 ; ...
                  500     6 ; ...
                  600    -7 ];

  % Simulate the test dataset

  y=sigma0*randn(n,1);
  y(simoutliers(:,1))=simoutliers(:,2);
  for k=1:size(simsteps,1)
    s=zeros(size(y));
    s(simsteps(k,1):end,1)=1;
    y=y+simsteps(k,2)*s;
  end

  % Set input parameters for function

  nwindow=11;
  varargin={ 'verbose', 2, 'simoutliers', simoutliers ,'simsteps', simsteps }

end

%% Input argument processing

% Check mandatory arguments

if ~isempty(mfilename) && nargin < 2, error('too few arguments (specify at least Y and NWINDOW)');end

% Default options

opt.verbose=0;
opt.crit=3.29;
opt.stepcrit=1.5;
opt.simoutliers=[];
opt.simsteps=[];
opt.simulation=false;

% Check options

for iopt=1:2:numel(varargin)
  optname=lower(varargin{iopt});
  if ~isfield(opt,optname), error(['Invalid option ' varargin{iopt} ]); end
  opt.(optname)=varargin{iopt+1};
  if any(strcmp(optname,{'simoutliers','simsteps'})), opt.simulation=true; end 
end

% Prepare data array

y=y(:);
n=size(y,1);

% Set window length (must be uneven) 

nwindow2=floor(nwindow/2);
nwindow=nwindow2*2+1;

%% Moving median computation with MAD estimate of standard deviation,

% The input dataset is replicated into a matrix ylagmatrix with nwindow columns,
% with in each column a shifted version of the dataset, padded by NaN's,
% so each row contains the data for the moving median computation.  

ylagmatrix=nan(n+nwindow-1,nwindow);
for k=1:nwindow
  ylagmatrix(k:n+k-1,k)=y;
end
ymm=nanmedian(ylagmatrix,2);
symm=1.4826*nanmedian( abs( ylagmatrix - repmat(ymm,[1 nwindow]) ), 2);
nymm=sum(~isnan(ylagmatrix),2);

% robust estimation of overall standard deviation (MAD)

sytot=1.4826*median(abs(y-ymm(nwindow2+1:end-nwindow2)));

% Outlier detection:
%
% For outlier detection we compute a matrix with residuals of the time
% series with the moving median. The matrix with residuals is reshaped, so 
% that rows are aligned on the datapoints.

yresmatrix=ylagmatrix - repmat(ymm,[1 nwindow]);
yresmatrix=reshape([nan(nwindow,1) ; yresmatrix(:)],[n+nwindow,nwindow]);
yresmatrix(1:nwindow,:)=[];

% Now we can do the outlier detection and count the number of moving
% windows in which an outlier is detected. Outliers typically have a score
% close to the window length. On the other hand, steps, result in scores
% that are less than half the window length.

% Compute correction factor f for the median residual. The standard deviation of
% the median is sigma_median= 1.253/sqrt(nwindow)*sigma. We use for nwindow
% half the window length and for sigma sytot. The moving MAD is not very
% usuable because of its inaccuracy, so we don't use it, and use sytot
% instead.

f=sqrt((nwindow2+1.253^2)/nwindow2);   

outliervote=sum( abs( yresmatrix  ) > f*sytot*opt.crit , 2);
idxoutlier=find(outliervote > nwindow2);
idxothervote=find(outliervote > 0 & outliervote <= nwindow2);

statoutlier=[  y(idxoutlier)-ymm(nwindow2+idxoutlier)  ...
               (y(idxoutlier)-ymm(nwindow2+idxoutlier))./(f*sytot) ... 
               outliervote(idxoutlier) ];

% Print

if opt.verbose >= 1
  fprintf('\nDetected outliers: %d   (sy=%.3f, f*sy=%.3f, limit=%.3f)\n\n', ...
      numel(idxoutlier),sytot,f*sytot,f*sytot*opt.crit)
  if numel(idxoutlier) > 0, fprintf('   idx      value normalized votes\n'); end
  for k=1:numel(idxoutlier)
    fprintf('%6d %10.3f %10.2f  %3d\n',idxoutlier(k), ...
            y(idxoutlier(k))-ymm(nwindow2+idxoutlier(k)), ...
            (y(idxoutlier(k))-ymm(nwindow2+idxoutlier(k)))./(f*sytot), ...
            outliervote(idxoutlier(k)));
  end
end

% Plot the series and median filtering result

if opt.verbose >= 2
  figure;
  plot(y,'.');
  hold on;
  plot(ymm(nwindow2+1:end-nwindow2),'g-','linewidth',2)
  plot(symm(nwindow2+1:end-nwindow2),'r:','linewidth',2)
  plot(idxoutlier,y(idxoutlier),'rx')
  text(idxoutlier,y(idxoutlier),num2str(outliervote(idxoutlier)),'VerticalAlignment','bottom');
  plot(idxothervote,y(idxothervote),'mo')
  text(idxothervote,y(idxothervote),num2str(outliervote(idxothervote)),'VerticalAlignment','top','Color','m');
  legendtext={'data' , 'median', 'stdev MAD', 'detected outliers' , 'outlier votes'};
  if opt.simulation
    plot(opt.simoutliers(:,1),y(opt.simoutliers(:,1)),'rd')
    limits=ylim();
    for k=1:size(opt.simsteps,1)
      line([opt.simsteps(k)  ; opt.simsteps(k) ],limits,'color','r','linestyle',':');
    end
    legendtext=[legendtext { 'simulated outliers' , 'simulated steps'  }];
  end
  legend(legendtext)
  title([ 'Moving median filter with outlier detection (window=' num2str(nwindow) ')' ] )
end

%% Moving average computation with estimate of standard deviation

% Moving average with outliers not removed

yma0=nanmean(ylagmatrix,2);     
syma0=nanstd(ylagmatrix,0,2);

% With outliers - detected during median filtering - removed

for k=1:nwindow
  ylagmatrix(idxoutlier+k-1,k)=nan;
end
yma=nanmean(ylagmatrix,2);     
syma=nanstd(ylagmatrix,0,2);
nyma=sum(~isnan(ylagmatrix),2);

% Plot the moving averages

if opt.verbose >= 2
  figure;
  plot(y,'.');
  hold on;
  plot(ymm(nwindow2+1:end-nwindow2),'g-','linewidth',2)
  plot(symm(nwindow2+1:end-nwindow2),'r:','linewidth',2)
  plot(yma0(nwindow2+1:end-nwindow2),'m-','linewidth',1)
  plot(syma0(nwindow2+1:end-nwindow2),'m:','linewidth',1)
  plot(yma(nwindow2+1:end-nwindow2),'k-','linewidth',2)
  plot(syma(nwindow2+1:end-nwindow2),'k:','linewidth',2)
  plot(idxoutlier,y(idxoutlier),'rx')
  text(idxoutlier,y(idxoutlier),num2str(outliervote(idxoutlier)),'VerticalAlignment','bottom');
  legendtext={'data' , 'median', 'stdev MAD', 'mean (w/ outliers)' , 'stdev (w/ outliers)' , 'mean (w/o outlier)' , 'stdev (w/o outlier)' , 'detected outliers'};
  if opt.simulation
    plot(opt.simoutliers(:,1),y(opt.simoutliers(:,1)),'rd')
    limits=ylim();
    for k=1:size(opt.simsteps,1)
      line([opt.simsteps(k)  ; opt.simsteps(k) ],limits,'color','r','linestyle',':');
    end
    legendtext=[legendtext { 'simulated outliers' , 'simulated steps'  }];
  end
  legend(legendtext);
  title([ 'Moving median and mean with error estimates (window=' num2str(nwindow) ')' ] )
end

%% Bail out if no step detection is to be done

if isnan(opt.stepcrit) || opt.stepcrit <= 0 
  ym=yma(nwindow2+1:end-nwindow2);
  sy=syma(nwindow2+1:end-nwindow2);
  ny=nyma(nwindow2+1:end-nwindow2);
  idxstep=[];
  steps=[];
  return
end

%% Step detection

% Compute discriminator for the step detection. The step itself is estimated
% as the the difference of the moving average/median before and after the step. 
% The discriminator is the normalized valued which can be used for a test.
%
% There are several posibilities:
%
% - use the moving median with overlapping windows, and scale with the MAD.
% - use the moving median, non-overlapping windows, and scale with standard 
%   deviation (after outlier removal)
% - use the moving mean (after outlier removal, non-overlapping windows,
%   and scale with standard deviation (after outlier removal)
% - use switching edge detector (SED)
%
% If we would use the median with non-overlapping windows then nwindow2
% datapoints around the step will be identified (block shape). Scaling with 
% the standard deviation will result in a triangular shaped function. However,
% if we scale with the MAD, then the block shape remains as the MAD is not
% very sensitive outliers from "the other side of the step". The is
% apparent from the following plot (commented out)
%
% figure
% for overlap2=0:nwindow2
%   ystep=ymm(nwindow-overlap2+1:end-overlap2)-ymm(1+overlap2:end-nwindow+overlap2);
%   ystepo=( ymm(nwindow-overlap2+1:end-overlap2) - ymm(1+overlap2:end-nwindow+overlap2) ) ./ ...
%          sqrt( symm(nwindow-overlap2+1:end-overlap2).^2 +symm(1+overlap2:end-nwindow+overlap2).^2 );
%   ystepo2=( ymm(nwindow-overlap2+1:end-overlap2) - ymm(1+overlap2:end-nwindow+overlap2) ) ./ ...
%          sqrt( syma(nwindow-overlap2+1:end-overlap2).^2 +syma(1+overlap2:end-nwindow+overlap2).^2 );
%   plot(2:n,ystep-overlap2*3,':',2:n,ystepo-overlap2*3,'--',2:n,ystepo2-overlap2*3,'-')
%   hold on
% end

% Moving median with overlapping windows, and scale with the MAD (just for visualization)

ystepmm2=[NaN ; ymm(nwindow-nwindow2+1:end-nwindow2)-ymm(1+nwindow2:end-nwindow+nwindow2) ];
dstepmm2=[NaN ; ( ymm(nwindow-nwindow2+1:end-nwindow2)-ymm(1+nwindow2:end-nwindow+nwindow2) ) ./ ...
         sqrt( symm(nwindow-nwindow2+1:end-nwindow2).^2 + symm(1+nwindow2:end-nwindow+nwindow2).^2 )]; 
dstepmm2a=[ NaN ; ( ymm(nwindow-nwindow2+1:end-nwindow2)-ymm(1+nwindow2:end-nwindow+nwindow2) ) ./ sqrt(2)*sytot ];

% Moving median, non-overlapping windows, and scale with standard deviation (after outlier removal)

ystepmm=[NaN ; ymm(nwindow+1:end) - ymm(1:end-nwindow) ];
dstepmm=[NaN ; ( ymm(nwindow+1:end) - ymm(1:end-nwindow) ) ./ ...
         sqrt( syma(nwindow+1:end).^2 + syma(1:end-nwindow).^2 ) ];

idxstepmm=find( abs(dstepmm) > opt.stepcrit & ...
                abs(dstepmm) > [abs(dstepmm(2:end)) ; 0 ]       & ...
                abs(dstepmm) > [0 ; abs(dstepmm(1:end-1)) ]     & ...
                abs(dstepmm) > [abs(dstepmm(3:end)) ; 0 ; 0 ]   & ...
                abs(dstepmm) > [0 ; 0 ; abs(dstepmm(1:end-2))]  );

% Moving mean (after outlier removal), non-overlapping windows, and scale with standard deviation (after outlier removal)

ystepma=[ NaN ; yma(nwindow+1:end)-yma(1:end-nwindow) ];
dstepma=[ NaN ; ( yma(nwindow+1:end) - yma(1:end-nwindow) ) ./ ...
         sqrt( syma(nwindow+1:end).^2 + syma(1:end-nwindow).^2 ) ];

idxstepma=find( abs(dstepma) > opt.stepcrit & ...
                abs(dstepma) > [abs(dstepma(2:end)) ; 0 ]       & ...
                abs(dstepma) > [0 ; abs(dstepma(1:end-1)) ]     & ...
                abs(dstepma) > [abs(dstepma(3:end)) ; 0 ; 0 ]   & ...
                abs(dstepma) > [0 ; 0 ; abs(dstepma(1:end-2))]  );

% Print

if opt.verbose >= 1
  fprintf('\nDetected steps moving median method: %d\n\n',numel(idxstepmm))
  if numel(idxstepmm) > 0, fprintf('   idx      value normalized\n'); end
  for k=1:numel(idxstepmm)
    fprintf('%6d %10.3f %10.2f\n',idxstepmm(k),ystepmm(idxstepmm(k)),dstepmm(idxstepmm(k)));
  end
  fprintf('\nDetected steps moving average method: %d\n\n',numel(idxstepma))
  if numel(idxstepma) > 0, fprintf('   idx      value normalized\n'); end
  for k=1:numel(idxstepma)
    fprintf('%6d %10.3f %10.2f\n',idxstepma(k),ystepma(idxstepma(k)),dstepma(idxstepma(k)));
  end
end

% Plot

if opt.verbose >= 2
  figure
  subplot(2,1,1)
  plot(1:n,ystepmm,'k-',1:n,ystepma,'b-',1:n,ystepmm2,'g-')
  limits=ylim();
  if opt.simulation
    for k=1:size(opt.simsteps,1)
      line([opt.simsteps(k)  ; opt.simsteps(k) ],limits,'color','r','linestyle',':');
    end
  end
  for k=1:numel(idxstepmm)
    line([idxstepmm(k)  ; idxstepmm(k) ],limits,'color','g','linestyle','--');
    %text(idxstepmm(k),limits(1),[ ' mmstep=',num2str(ystepmm(idxstepmm(k))) ],'verticalalignment','bottom')
  end
  for k=1:numel(idxstepma)
    line([idxstepma(k)  ; idxstepma(k) ],limits,'color','r','linestyle','--');
    %text(idxstepma(k),limits(2),[ ' mastep=' num2str(ystepma(idxstepma(k))) ],'verticalalignment','top')
  end
  title([ 'Estimated step size (window=' num2str(nwindow) ')' ] )
  legend({'moving median' , 'moving average', 'overlapped moving median' })
  ylabel('[physical units]')
  subplot(2,1,2)
  plot(1:n,dstepmm,'k-',1:n,dstepma,'b-',1:n,dstepmm2,'g-',1:n,dstepmm2a,'r:')
  limits=ylim();
  if opt.simulation
    for k=1:size(opt.simsteps,1)
      line([opt.simsteps(k)  ; opt.simsteps(k) ],limits,'color','r','linestyle',':');
    end
  end
  for k=1:numel(idxstepmm)
    line([idxstepmm(k)  ; idxstepmm(k) ],limits,'color','g','linestyle','--');
  end
  for k=1:numel(idxstepma)
    line([idxstepma(k)  ; idxstepma(k) ],limits,'color','r','linestyle','--');
  end
  limits=xlim();
  line(limits,[opt.stepcrit ; opt.stepcrit],'color','r','linestyle',':');
  line(limits,[-opt.stepcrit ; -opt.stepcrit],'color','r','linestyle',':');
  title([ 'Step discriminator function (window=' num2str(nwindow) ')' ] )
  ylabel('[-]')
  legend({'moving median' , 'moving average', 'overlapped moving median' })
end


%% Switching Edge Detector (SED)
%
% Smith (1998) introduced the switching edge detector (SED) for detecting 
% steps in time series using moving averages. The principle is the same,
% two windows, before and after a supposed step, and the test quantity
% normalized by the emperical standard deviation. What is different, is
% that two weighting factors are introduced, gm (before) and gp (after),
% which implement a steep switching function.
%
% The weighting functions are computed below

r=50;
r2=r*2;
gm= syma(1:end-nwindow).^r2  ./ ( syma(nwindow+1:end).^r2 + syma(1:end-nwindow).^r2 );
gp= syma(nwindow+1:end).^r2  ./ ( syma(nwindow+1:end).^r2 + syma(1:end-nwindow).^r2 );

% This results in the "switching" moving average ysed and sed detector function

ysed=[NaN ; gm.*yma(nwindow+1:end)  + gp.*yma(1:end-nwindow)];

dsed=[ NaN ; ( yma(nwindow+1:end) - yma(1:end-nwindow) ) ./ ...
         sqrt( gp.*syma(nwindow+1:end).^2 + gm.*syma(1:end-nwindow).^2 ) ];

idxstepsed=find( abs(dsed) > opt.stepcrit & ...
                 abs(dsed) > [abs(dsed(2:end)) ; 0 ]   & ...
                 abs(dsed) > [0 ; abs(dsed(1:end-1)) ]);

% Print

if opt.verbose >= 1
  fprintf('\nDetected steps switching edge detector method (SED): %d\n\n',numel(idxstepsed))
  if numel(idxstepsed) > 0, fprintf('   idx      value normalized\n'); end
  for k=1:numel(idxstepsed)
    fprintf('%6d %10.3f %10.2f\n',idxstepsed(k),ystepma(idxstepsed(k)),dsed(idxstepsed(k)));
  end
end

% Plot

if opt.verbose >= 2
  figure
  subplot(2,1,1)
  plot(1:n,yma(nwindow2+1:end-nwindow2),'k-',1:n,ymm(nwindow2+1:end-nwindow2),'g-',1:n,ysed,'b-')
  limits=ylim();
  if opt.simulation
    for k=1:size(opt.simsteps,1)
      line([opt.simsteps(k)  ; opt.simsteps(k) ],limits,'color','r','linestyle',':');
    end
  end
  for k=1:numel(idxstepsed)
    line([idxstepsed(k)  ; idxstepsed(k) ],limits,'color','g','linestyle','--');
    %text(idxstepmm(k),limits(1),[ ' mmstep=',num2str(ystepmm(idxstepmm(k))) ],'verticalalignment','bottom')
  end
  title([ 'Moving mean, median and SED estimate (window=' num2str(nwindow) ')' ] )
  legend({'moving mean', 'moving median' , 'SED' })
  ylabel('[physical units]')
  subplot(2,1,2)  
  plot(1:n,dstepma,'k-',1:n,dstepmm,'g-',1:n,dsed,'b-')
  limits=ylim();
 if opt.simulation
    for k=1:size(opt.simsteps,1)
      line([opt.simsteps(k)  ; opt.simsteps(k) ],limits,'color','r','linestyle',':');
    end
  end
  for k=1:numel(idxstepsed)
    line([idxstepsed(k)  ; idxstepsed(k) ],limits,'color','g','linestyle','--');
    %text(idxstepmm(k),limits(1),[ ' mmstep=',num2str(ystepmm(idxstepmm(k))) ],'verticalalignment','bottom')
  end
  limits=xlim();
  line(limits,[opt.stepcrit ; opt.stepcrit],'color','r','linestyle',':');
  line(limits,[-opt.stepcrit ; -opt.stepcrit],'color','r','linestyle',':');
  title([ 'Step discriminator function (window=' num2str(nwindow) ')' ] )
  legend({'moving mean', 'moving median' , 'SED' })
  ylabel('[-]')
end


%% Correct moving average for steps

% Select final method for step detection

idxstep=idxstepsed;
steps=ystepma(idxstepsed);

% Apply steps to ylagmatrix and recompute the moving average

for l=1:numel(idxstep)
  ll=idxstep(l);
  lll=ll+nwindow2-1;
  for k=1:nwindow2
    ylagmatrix(ll+k-1:lll,k)=ylagmatrix(ll+k-1:lll,k)-steps(l);
  end
  for k=nwindow2+2:nwindow
    ylagmatrix(lll+1:ll+k-2,k)=ylagmatrix(lll+1:ll+k-2,k)+steps(l);
  end
end
yma2=nanmean(ylagmatrix,2);     
syma2=nanstd(ylagmatrix,0,2);
nyma2=sum(~isnan(ylagmatrix),2);

% Plot

if opt.verbose >= 2
  figure;
  plot(y,'.');
  hold on;
  plot(ymm(nwindow2+1:end-nwindow2),'g-','linewidth',2)
  plot(symm(nwindow2+1:end-nwindow2),'r:','linewidth',2)
  plot(yma0(nwindow2+1:end-nwindow2),'m-','linewidth',1)
  plot(syma0(nwindow2+1:end-nwindow2),'m:','linewidth',1)
  plot(yma(nwindow2+1:end-nwindow2),'k-','linewidth',2)
  plot(syma(nwindow2+1:end-nwindow2),'k:','linewidth',2)
  plot(yma2(nwindow2+1:end-nwindow2),'b-','linewidth',2)
  plot(syma2(nwindow2+1:end-nwindow2),'b:','linewidth',2)
  plot(idxoutlier,y(idxoutlier),'rx')
  text(idxoutlier,y(idxoutlier),num2str(outliervote(idxoutlier)),'VerticalAlignment','bottom');
  legendtext={'data' , 'median', 'stdev MAD', 'mean (w/ outliers)' , 'stdev (w/ outliers)' , 'mean (w/o outlier)' , 'stdev (w/o outlier)' , 'robust mean' , 'robust stdev'  };
  limits=ylim();
  for k=1:numel(idxstep)
    line([idxstep(k)  ; idxstep(k) ],limits,'color','g','linestyle','--');
    % text(idxstep(k),limits(2),[ ' step=' num2str(step(k)) ],'verticalalignment','top')
  end
  if opt.simulation
    plot(opt.simoutliers(:,1),y(opt.simoutliers(:,1)),'rd')
    for k=1:size(opt.simsteps,1)
      line([opt.simsteps(k)  ; opt.simsteps(k) ],limits,'color','r','linestyle',':');
    end
    %legendtext=[legendtext { 'simulated outliers' , 'simulated steps'  }];
  end
  legend(legendtext);
  title([ 'Moving median and mean with error estimates (window=' num2str(nwindow) ')' ] )
end

%%

% reduce ym and ys to length of the input time series (with nxtr points at
% the start and end of the series)

ym=yma2(nwindow2+1:end-nwindow2);
sy=syma2(nwindow2+1:end-nwindow2);
ny=nyma2(nwindow2+1:end-nwindow2);

end