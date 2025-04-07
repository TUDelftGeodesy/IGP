function varargout=stmresiduals(pntName, pntCrd, epochIds, dslocal, lsqstat1, opt, outputId)
%stmresiduals   Residual Testing for Integrated processing..
%   STMRESIDUALS(pntName,pntCrd,epochIds,dslocal,lsqstat1,opt,outputId) tests
%   the normalized least squares residuals. The following test statistics are
%   computed
%   - w-test statistic for individual observations
%   - overall model test statistic for each point/epoch set (of displacements)
%   - overall model test statistic per point
%   - overall model test statistic per epoch
%   - overall model test statistic per dataset
%   - global overall model test statistic
%   The observations and points with the largest test statistics are 
%   printed, including all epoch and dataset test statistics. 
%
%   Input is the cell array dslocal that is computed by stmintegrate. Other
%   inputs are the pntNames, point coordinates (pntNeu) and epochIds. Lsqstat1
%   contains the parameter and observation counts. Options from
%   stmintegrate are passed in opt, with 
%
%      opt.alpha0=0.001;          % probability of false alarm
%
%      opt.critWtest=nan;         % criterion for printing test statistics,
%      opt.critDisplOmt=nan;      % if nan, then the criterion is computed
%      opt.critPointOmt=nan;      % from the probability of false alarm
%      opt.critEpochOmt=nan;      % opt.alpha and the redundancy
%
%      opt.maxWtest=50;           % maximum of test statistics to print
%      opt.maxDisplOmt=20;
%      opt.maxPointOmt=10;
%
%   If opt.doplots > 0 then histograms of the residuals and point statistics
%   are plotted, using outputId as title. 
%
%   Examples:
%
%      stmresiduals(pntName, pntCrd, epochIds, dslocal, lsqstat1, opt, datasetId)
%
%      [omtPoints,dofPoints,omtEpochs,dofEpochs,omtDatasets,dofDatasets]= ...
%          stmresiduals(pntName, pntCrd, epochIds, dslocal, lsqstat1, opt, datasetId)
%
%      outputfilename='groningen_min1peryear_U22combined.mat';
%      [filepath,outputId] = fileparts(outputfilename);
%      load(fullfile(filepath,[ outputId '_res.mat'])); % read residuals from mat file
%      opt.doplots=1; % turn plotting on
%      stmresiduals(pntName, pntCrd, epochIds, dslocal, lsqstat1, opt, outputId)
%
%   The output omtPoints and omtEpochs have already the pntMask and
%   epochMask applied...
%
%   See also STMINTEGRATE and STMRESPLOT.
%
%  (c) Hans van der Marel, Delft University of Technology, 2021.

% Created:  13 February 2021 by Hans van der Marel
% Modified: 02 March 2021 by Hans van der Marel
%              - initial release              
%           02 June 2022 by Hans van der Marel
%              - replaced pntId by pntName, adjusted formats  
%              - use pntCrd instead of pntNeu, if opt.doplots_mapcrd=true,
%                then pntCrd is expected to contain map coordinates
%              - plotlabels as variable
%           26 October 2022 by Hans van der Marel
%              - increased length of pntName in print statements (28 -> 38)
%           13 April 2023 by Hans van der Marel
%              - fixed bug in overlapping residual plot plotting residuals
%                instead of normalized residuals
%           21 November 2023 by Hans van der Marel
%              - added background and defoborder to map plots

%% Initialize 

% Number of datasets, points and epochs

numDatasets=numel(dslocal);   
numPoints=numel(pntName);
numEpochs=numel(epochIds);

% Compute point and epoch mask

pntMask=false(size(pntName));
epochMask=false(size(epochIds));
for k=1:numDatasets  
    for l=1:dslocal{k}.numObsTypes
       pntMask(dslocal{k}.lsq{l}.idxPnt)=true;
       epochMask(dslocal{k}.lsq{l}.idxEpo)=true;
    end
end

% Convert masks to index arrays

idxPointMask=1:numel(pntMask);
idxPointMask=idxPointMask(pntMask);

idxEpochMask=1:numel(epochMask);
idxEpochMask=idxEpochMask(epochMask);

   
%% Print individual outliers (normalized residuals exceeding criterion)

v=1;
if isnan(opt.critWtest)
   critWtest = sqrt(chi2inv(1-opt.alpha0,v)/v);
else
   critWtest = opt.critWtest;
end

outlierWtest=[];
for k=1:numDatasets  
    for l=1:dslocal{k}.numObsTypes
       idxmask=dslocal{k}.lsq{l}.idxmask;
       en=dslocal{k}.lsq{l}.en;
       e=dslocal{k}.lsq{l}.e;
       ien=find(abs(en) > critWtest);
       if numel(ien) > 0
          outlierEpochs=floor(idxmask(ien)./numPoints);
          outlierPoints=idxmask(ien)-outlierEpochs.*numPoints;      
          outlierWtest = [ outlierWtest ; ...
              repmat(k,size(ien)) repmat(l,size(ien)) outlierPoints outlierEpochs en(ien) e(ien) ];
       end
    end
end
if size(outlierWtest,1) > 0
   [~,ien]=sort(abs(outlierWtest(:,5)),'descend');
   fprintf('\nOutliers normalized residuals (critical value %.3f, shown max %d rejections):\n\n',critWtest,opt.maxWtest)
   fprintf('pntName                                epochId      datasetId                              obsType  normRes    Res\n')
   for k=1:min(size(outlierWtest,1),opt.maxWtest)
       kk=ien(k);
       fprintf('%-38s %-12s %-39s %-8s %6.2f %6.2f\n',pntName{outlierWtest(kk,3)},epochIds{outlierWtest(kk,4)},dslocal{outlierWtest(kk,1)}.datasetId,dslocal{outlierWtest(kk,1)}.obsTypes{outlierWtest(kk,2)},outlierWtest(kk,5),outlierWtest(kk,6))
   end
end

%% Compute OMT per point and epoch 
%
% Compute OMT (reduced Chi-square) like statistics
%
% - reduced Chi-square, i.e. divided by degree of freedom, i.e. expected value is one  (also Fisher-Snedecor F(dof,inf,0) distribution)
% - take the square root
% - the (absolute value of) W-test is the one dimensional version of this test 
%
% OMT is not a very good name for these, it is actual the square root of 
% the "mean squared weighted deviation (MSWD)", or sMSWD

% Collect OMT data in point/epoch matrix and compute degree of freedom

omtDispl=zeros(numPoints,numEpochs);
dofDispl=zeros(numPoints,numEpochs);

for k=1:numDatasets  
    for l=1:dslocal{k}.numObsTypes
       en=dslocal{k}.lsq{l}.en;
       rn=dslocal{k}.lsq{l}.rn;
       idxmask=dslocal{k}.lsq{l}.idxmask;
       omtDispl(idxmask)=omtDispl(idxmask)+en.^2;
       dofDispl(idxmask)=dofDispl(idxmask)+rn;
    end
end
dofDispl2=dofDispl(pntMask,epochMask);

omtDispl2=sqrt(omtDispl(pntMask,epochMask)./dofDispl2);
omtDispl2(dofDispl2<=0)=0;

v=median(dofDispl2(:));
if isnan(opt.critDisplOmt)
   critDisplOmt = sqrt(chi2inv(1-opt.alpha0,v)/v);
else
   critDisplOmt= opt.critDisplOmt;
end

[~,ien]=sort(omtDispl2(:),'descend');
fprintf('\nPoint/Epoch OMT (critical value %.3f, dof=%.1f, shown max %d rejections):\n\n',critDisplOmt,v,opt.maxDisplOmt);
fprintf('pntName                                epochId            dof       omt\n')
for k=1:opt.maxDisplOmt
  kk=ien(k);
  [kp,ke]= ind2sub(size(omtDispl2),kk);
  kp=idxPointMask(kp);
  ke=idxEpochMask(ke);
  if omtDispl2(kk) < critDisplOmt, break; end
  fprintf('%-38s %-12s %9.1f %9.2f\n',pntName{kp},epochIds{ke},dofDispl2(kk),omtDispl2(kk))
end


%% Print point OMT

dofPoints=sum(dofDispl2,2);
omtPoints=sqrt(sum(omtDispl(pntMask,epochMask),2)./dofPoints);
omtPoints(dofPoints<=0)=0;
dofPoints(dofPoints<=0)=0;

v=median(dofPoints(:));
if isnan(opt.critPointOmt)
   critPointOmt = sqrt(chi2inv(1-opt.alpha0,v)/v);
else
   critPointOmt= opt.critPointOmt;
end

[~,ien]=sort(omtPoints,'descend');
fprintf('\nPoint OMT (critical value %.3f, dof %.1f, shown max %d rejections):\n\n',critPointOmt,v,opt.maxPointOmt);
fprintf('pntName                                      dof  omtPoint\n')
for k=1:opt.maxPointOmt
  kk=ien(k);
  kp=idxPointMask(kk);
  if omtPoints(kk) < critPointOmt, break; end
  fprintf('%-38s %9.1f %9.2f\n',pntName{kp},dofPoints(kk),omtPoints(kk))
end

%% Plot point OMT (optional)

if opt.doplots > 0

    % Check if datatips are supported

    opt.datatips=true;
    if opt.datatips && exist('dataTipTextRow') ~= 2
        fprintf('Warning: Datatips are not supported (by this Matlab version)')
        opt.datatips=false;
    end
    
    % Map labels and aspect ratio
    
    if isfield(opt,'doplots_mapcrd') && opt.doplots_mapcrd
       xylabel = { 'East [km]', 'North [km]' };
       xyaspect = [1 1 1];
    else
       xylabel = { 'Longitude [deg]', 'Latitude [deg]' };
       xyaspect = [ 1/cosd(mean(pntCrd(:,1))) 1 1];
    end

    % Get background and deformation border

    %    'background'         Map background (Default 'igp-background.mat')
    %    'unstablearea'       Polygon defining unstable area (Default 'igp-unstablearea.mat')

    if ~isfield(opt,'background')
        opt.background='igp-background.mat';     % map background
    end
    if ~isfield(opt,'unstablearea')
        opt.defoborder='igp-unstablearea.mat'; % polygon defining unstable area
    end
    if ~isfield(opt,'defoborder')
        opt.defoborder='igp-unstablearea.mat'; % polygon defining unstable area
    end

    if isnumeric(opt.background)
        igpbackground=opt.background;
    elseif exist(opt.background,'file')
        igpbackground=load(opt.background);
        if isstruct(igpbackground)
            igpbackground = [ igpbackground.Lat(:) igpbackground.Lon(:) ];
        end
    else
        igpbackground=[];
    end
    
    if isnumeric(opt.defoborder)
        defoborder=opt.defoborder;
    elseif exist(opt.defoborder,'file')
        defoborder=load(opt.defoborder);
        if isstruct(defoborder)
            defoborder = [ defoborder.Lat(:) defoborder.Lon(:) ];
        end
    else
        defoborder=[];
    end

    
    % Pseudo color plot with point/epoch omt 

    figure;
    imagesc(sqrt(omtDispl2));
    ylabel('Points')
    xlabel('Epochs')
    gcb=colorbar;
    ylabel(gcb,'\surd omt','interpreter','tex');
    title('OMT per point and epoch')

    % Plot point omt in ascending order 

    figure;
    [omtPointsSorted,kk]=sort(omtPoints);
    hs0=plot(omtPointsSorted,'.');
    ylabel('\surd omt','interpreter','tex');
    xlabel('Points')
    title('OMT per point (in ascending order)')
    
    if opt.datatips
        htip2 = dataTipTextRow('dof',dofPoints(kk),'%.1f');
        htip3 = dataTipTextRow('id',pntName(idxPointMask(kk)));
        hs0.DataTipTemplate.DataTipRows(1).Label = 'seqCount';
        hs0.DataTipTemplate.DataTipRows(2).Label = 'sOMT';
        hs0.DataTipTemplate.DataTipRows(end+1) = htip2;
        hs0.DataTipTemplate.DataTipRows(end+1) = htip3;
    end

    % Map with point omt and dof 

    figure;
    hs1=scatter(pntCrd(pntMask,2),pntCrd(pntMask,1),1+omtPoints*16/max(omtPoints),omtPoints);
    xlabel(xylabel{1});
    ylabel(xylabel{2});
    gcb=colorbar;
    ylabel(gcb,'\surd omt','interpreter','tex');
    title('Point OMT')
    hold on
    xl=xlim();yl=ylim(); 
    if ~isempty(igpbackground)
        plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
    end
    if ~isempty(defoborder)
        plot(defoborder(:,2),defoborder(:,1),':m')
    end
    xlim(xl);ylim(yl);
    daspect(xyaspect);
    
    if opt.datatips
        htip1 = dataTipTextRow('OMT',omtPoints,'%.2f');
        htip2 = dataTipTextRow('dof',dofPoints,'%.1f');
        htip3 = dataTipTextRow('id',pntName(idxPointMask(:)));
        hs1.DataTipTemplate.DataTipRows(1).Label = xylabel{1};
        hs1.DataTipTemplate.DataTipRows(2).Label = xylabel{2};
        hs1.DataTipTemplate.DataTipRows(end+1) = htip1;
        hs1.DataTipTemplate.DataTipRows(end+1) = htip2;
        hs1.DataTipTemplate.DataTipRows(end+1) = htip3;
    end
    
    figure;
    hs2=scatter(pntCrd(pntMask,2),pntCrd(pntMask,1),1+dofPoints*16/max(dofPoints),dofPoints);
    xlabel(xylabel{1});
    ylabel(xylabel{2});
    gcb=colorbar;
    ylabel(gcb,'dof');
    title('Degree of freedom Point OMT')   
    hold on
    xl=xlim();yl=ylim(); 
    if ~isempty(igpbackground)
        plot(igpbackground(:,2),igpbackground(:,1),'Color',[.5 .5 .5])
    end
    if ~isempty(defoborder)
        plot(defoborder(:,2),defoborder(:,1),':m')
    end
    xlim(xl);ylim(yl);
    daspect(xyaspect);
    
    if opt.datatips
        hs2.DataTipTemplate.DataTipRows(1).Label = xylabel{1};
        hs2.DataTipTemplate.DataTipRows(2).Label = xylabel{2};
        hs2.DataTipTemplate.DataTipRows(end+1) = htip1;
        hs2.DataTipTemplate.DataTipRows(end+1) = htip2;
        hs2.DataTipTemplate.DataTipRows(end+1) = htip3;
    end

end

%% Print epoch OMT

dofEpochs=sum(dofDispl2,1);
omtEpochs=sqrt(sum(omtDispl(pntMask,epochMask),1)./dofEpochs);
omtEpochs(dofEpochs<=0)=0;

v=median(dofEpochs(:));
if isnan(opt.critEpochOmt)
   critEpochOmt = sqrt(chi2inv(1-opt.alpha0,v)/v);
else
   critEpochOmt = opt.critEpochOmt;
end

[~,ien]=sort(omtEpochs,'descend');
fprintf('\nEpoch OMT (critical value %.3f, dof %.1f):\n\n',critEpochOmt,v);
fprintf('epochId                dof  omtEpoch\n')
rej='suspect';
for k=1:numel(omtEpochs)
  kk=ien(k);
  ke=idxEpochMask(kk);
  if omtEpochs(kk) < critEpochOmt, rej=''; end
  fprintf('%-16s %9.1f %9.2f  %s\n',epochIds{ke},dofEpochs(kk),omtEpochs(kk),rej)
end

%% Print dataset OMT

omtDatasets=zeros(numDatasets,1);
dofDatasets=zeros(numDatasets,1);

fprintf('\nDataset OMT:\n\n');
fprintf('datasetId                              obsType         dof       omt\n')
for k=1:numDatasets  
    for l=1:dslocal{k}.numObsTypes
       en=dslocal{k}.lsq{l}.en;
       rn=dslocal{k}.lsq{l}.rn;
       omtDs=sum(en.^2);
       %npar2k=dslocal{k}.lsqsave{l}.np+dslocal{k}.lsqsave{l}.ne-1;
       %dofDs=numel(en)-npar2k-npar2k/npar2*npar1;
       dofDs=sum(rn);
       omtDatasets(k)=omtDatasets(k)+omtDs;
       dofDatasets(k)=dofDatasets(k)+dofDs;
       fprintf('%-39s %-8s %9.1f %9.2f\n',dslocal{k}.datasetId,dslocal{k}.obsTypes{l},dofDs,sqrt(omtDs/dofDs))
    end
    if dslocal{k}.numObsTypes > 1
       fprintf('%-39s %-8s %9.1f %9.2f\n','','all',dofDatasets(k),sqrt(omtDatasets(k)/dofDatasets(k)))
    end
end
fprintf('                                                 ---------  --------\n')
fprintf('All                                              %9.1f %9.2f     (%d, %.2f)\n\n',sum(dofDatasets),sqrt(sum(omtDatasets)/sum(dofDatasets)),lsqstat1.dof,sqrt(lsqstat1.omt))


%% Histogram of normalized residuals (overlapping)

if opt.doplots > 0

    % Full range
    
    figure
    subplot(3,1,[1,2])
    for k=1:numDatasets
        for l=1:dslocal{k}.numObsTypes
            histogram(dslocal{k}.lsq{l}.en,'normalization','prob','displayname',[dslocal{k}.datasetId ' (' dslocal{k}.obsTypes{l} ')'])
            hold on
        end
    end
    legend('interpreter','none')
    ylabel('Probability')
    hold off

    % Tails

    subplot(3,1,3)
    for k=1:numDatasets
        for l=1:dslocal{k}.numObsTypes
            histogram(dslocal{k}.lsq{l}.en,'normalization','count','displayname',[dslocal{k}.datasetId ' (' dslocal{k}.obsTypes{l} ')'])
            hold on
        end
    end
    %legend('interpreter','none')
    ylabel('Counts')
    ylim([0 10])
    hold off
    
    sgtitle(['Normalized residuals ' outputId ],'interpreter','none')

end

%% Histogram of normalized residuals (subplots per dataset)

if opt.doplots > 0

    figure
    numObsGroups=sum(cellfun(@(x) x.numObsTypes,dslocal));
    n=min(4,numObsGroups);
    m=ceil(numObsGroups/n);
    ll=0;
    for k=1:numDatasets  
        for l=1:dslocal{k}.numObsTypes
           ll=ll+1;
           subplot(m,n,ll)
           histogram(dslocal{k}.lsq{l}.en,'normalization','prob')
           ylabel([dslocal{k}.datasetId ' (' dslocal{k}.obsTypes{l} ')'],'interpreter','none')
           text(.05,.95,{ sprintf('n=%d',numel(dslocal{k}.lsq{l}.en));  sprintf('\\sigma=%.2f',std(dslocal{k}.lsq{l}.en)) },'units','normalized')
        end
    end
    sgtitle(['Normalized residuals ' outputId ],'interpreter','none')

end

%% Histogram of residuals (overlapping)

if opt.doplots > 0

    % Full range

    figure
    subplot(3,1,[1,2])
    for k=1:numDatasets  
        for l=1:dslocal{k}.numObsTypes
           e=dslocal{k}.lsq{l}.e;
           histogram(e,'normalization','prob','displayname',[dslocal{k}.datasetId ' (' dslocal{k}.obsTypes{l} ')'])
           hold on
        end
    end
    legend('interpreter','none')
    ylabel('Probability')
    hold off

    % Tails
    
    subplot(3,1,3)
    for k=1:numDatasets  
        for l=1:dslocal{k}.numObsTypes
           e=dslocal{k}.lsq{l}.e;
           histogram(e,'normalization','count','displayname',[dslocal{k}.datasetId ' (' dslocal{k}.obsTypes{l} ')'])
           hold on
        end
    end
    %legend('interpreter','none')
    ylabel('Counts')
    xlabel('mm')
    ylim([0 10])
    hold off

    sgtitle(['Residuals ' outputId ],'interpreter','none')

end

%% Histogram of residuals (subplots per dataset)

if opt.doplots > 0

    figure
    numObsGroups=sum(cellfun(@(x) x.numObsTypes,dslocal));
    n=min(4,numObsGroups);
    m=ceil(numObsGroups/n);
    ll=0;
    for k=1:numDatasets  
        for l=1:dslocal{k}.numObsTypes
           ll=ll+1;
           subplot(m,n,ll)
           e=dslocal{k}.lsq{l}.e;
           histogram(e,'normalization','prob')
           ylabel([dslocal{k}.datasetId ' (' dslocal{k}.obsTypes{l} ')'],'interpreter','none')
           xlabel('mm')
           text(.05,.95,{ sprintf('n=%d',numel(e));  sprintf('\\sigma=%.2f mm',std(e)) },'units','normalized')
        end
    end
    sgtitle(['Residuals ' outputId ],'interpreter','none')

end

if nargout > 0
   varargout = { omtPoints dofPoints omtEpochs dofEpochs omtDatasets dofDatasets};
end

end