function hf=stmplotseries(st,varargin)
%STMPLOTSERIES Plot timescomponents for multiple stations.
%   STMPLOTSERIES(ST) plots time series for stations in the space time matrix ST. 
%   ST can be a space time matrix structure, or the name of a space time matrix
%   dataset.
%
%   STMPLOTSERIES(ST,OPTION,VALUE,...) accepts multiple OPTION/VALUE pairs, valid 
%   option/value pairs are
%
%     item          string       Character string to select either 'obsData' or 'auxData' to plot (default 'obsData')
%     types         cell array   Cell array with selected obsTypes/auxType to plot (default all in 'obsData' or 'auxData')
%     events        cell array   Events to plot, default ({'ANT' 'DPL' 'STEPS'}), i.e no 'REC' 
%
%     sort          string       name of the point attribute to sort on (default none == order in space time matrix),
%                                must be a valid fieldname of the pntAttrib with subindexing, useful examples are 
%                                "vel(:,3)", "omt(:,3)", "rms(:,3)" and "amplitude(:,1,3)", where each time is
%                                sorted on the third (e.g. Up) component.
%
%     subplot       logical      If true use subplots (default), if false create individual plots, the
%                                number of subplots is equal to the number of observation types, or maxSubPlots, 
%                                or if one observation type is selected it depends on the number of stations and
%                                number of traces per subplot
%
%     maxPlots      integer      maximum number of plots (figures) (default 1) - cannot be changed at the moment
%     maxSubPlots   integer      maximum number of subplots per figure (default 4)
%     maxTraces     integer      maximum number of traces per plot/subplot (default 60)
%
%     spacing       scalar | 1x3 spacing between traces , default is 10 ([10 10 10])
%     rowMargin     scalar | 1x2 add rowmargin extra (blank) rows at top and bottom of plot (default auto -1),
%                                [bottom top], or scalar, if negative does auto scaling
%     position      1x4 array    figure Position [ left bottom width height] in pixels ([50 50 560 870]) 
%                                for the leftmost plot (or subplot)
%     innerPosition 1x4 array    axis InterPosition [ left bottom width height] in normalized units ([ 0.10 0.04 0.87 0.93]) 
%
%     scale         scalar       scale factor for data (default 1)
%     marker        logical      plot markers (default false) 
%
%     schranking    string       S-transform method (default 'none')
%                                  'none'    - none
%                                  'minnorm' - minumum norm solution
%                                  'sdnorm'  - minimum norm of single differences in time
%                                  'ssdnorm' - minimum norm of symmetric single differences  
%                                  'min-velocities' - S-transform based on velocity estimation
%     refPoints     cell array   Optional reference points for S-transform
%     unstablearea  string       Polygon defining unstable area (default 'igp-unstablearea.mat') 
%                                for S-transform 
%
%     saveplt       string        Directory for pdf plots (default [] is none)
%
%   HF=STMPLOTSERIES(...) returns the plot handles for the figures.
%
%   Examples:
%
%      stmplotseries(st)
%      stmplotseries(st,'item','aux')
%
%      stmplotseries(st,'sort','rms(:,3)','rowMargin',[2,4])
%      stmplotseries(st,'sort','vel(:,3)','rowMargin',[2,5])
%
%      stmplotseries(st,'types','Up','maxTraces',30,'rowMargin',[1 4])
%      stmplotseries(st,'types','Up','maxTraces',30,'rowMargin',[1 4],'sort','vel(:,3)')
%      stmplotseries(st,'types','Up','maxTraces',30,'rowMargin',[1 4],'sort','omt(:,3)')
%
%   See also STMDISP and STMDIFF.
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016,2021.

% Created:   5 April 2015 by Hans van der Marel
% Modified: 29 August 2016 by Hans van der Marel
%              - added to tseries toolbox
%            5 March 2021 by Hans van der Marel
%              - removed read file support
%              - changed the input arguments, alowing option/value pairs
%              - added options to change the layout
%              - renamed tseriesplotcomponents to tseriesplotcomp to avoid
%                backport issues
%           15 March 2021 by Hans van der Marel
%              - port of SERIESPLOTCOMP to STMPLOTSERIES
%           26 September 2023 by Hans van der Marel
%              - added options to plot markers and scale the data to mm
%              - fix for empty data rows
%            6 May 2024 by Hans van der Marel
%              - make sure epochs are sorted on time
%           16 May 2024 by Hans van der Marel
%              - S-transform (schranking) added
%           11 June 2024 by Hans van der Marel
%              - Added min-velocities method to S-transform
%           17 June 2024 by Hans van der Marel
%              - added option to save plots to pdf

% Default options

opt.item='obsData';         % select item to plot (default 'obsData')
opt.types='';               % Cell array with selected obsTypes/auxType to plot (default all in opt.item)
opt.events={'ANT' 'DPL' 'STEPS'};  % cell array with events to plot 

opt.spacing=10;                % spacing between traces
opt.rowMargin=-1;              % add rowmargin extra (blank) rows at top and bottom of plot (-1)
opt.subplot=true;              % Use subplots (default), if false create individual plots
opt.position=[50 50 560 870];   % figure Position [ left bottom width height] in pixels for the leftmost plot (or subplot)
%opt.innerPosition=[ 0.13 0.11 0.775 0.815];   % axis InnerPostion [ left bottom width height] in normilized units 
opt.innerPosition=[ 0.10 0.04 0.87 0.93];   % axis InnerPosition [ left bottom width height] in normalized units 

opt.sort='';                   % point attribute to sort on

opt.maxPlots = 1;              % maximum number of plots
opt.maxSubPlots = 4;           % maximum number of subplots
opt.maxTraces = 60;            % maximum number of traces per plot

opt.scale=1;                   % scale factor
opt.marker=false;              % if true, plot markers 
opt.xfmargin=4;                % excess margin in parts of year (TBC)
opt.FontSize=10;

opt.schranking='none';         % S-transformation method
opt.refPoints={};              % reference points for S-transform
opt.unstablearea='igp-unstablearea.mat'; % polygon defining unstable area

opt.saveplt=[];                % directory for pdf plots (default [] is none)

% Check arguments

if nargin < 1, error('too few arguments (specify at least stationames and item');end
if ischar(st)
   st=stmread(st);
elseif ~isstruct(st)
   error('The first argument must be a character string with the filename or space time matrix structure.')
end

% Check the option/value pairs

for k=1:2:numel(varargin)-1
   if isfield(opt,varargin{k})
      opt.(varargin{k})=varargin{k+1};
   else
      error(['Illegal option ' varargin{k} ])
   end
end

if isscalar(opt.spacing)
  spacing=[opt.spacing opt.spacing opt.spacing];
else
  spacing=opt.spacing;
end
if isscalar(opt.rowMargin)
  rowmargin=[opt.rowMargin opt.rowMargin];
else
  rowmargin=opt.rowMargin;
end

% Select items to plot

datasetId = st.datasetId;
pntName=st.pntName;
pntAttrib=st.pntAttrib;
year=st.epochDyear;

[year,iepoch]=sort(year);  % ensure year is sorted ...

item=opt.item;
if strncmpi(item,'obsData',3)
  obsTypes= st.obsTypes; 
  %obsData=st.obsData;
  obsData=st.obsData(:,iepoch,:);      % adjust to sorted epochs
elseif strncmpi(item,'auxData',3)
  fprintf('Plot auxiliary data types instead of observation types\n')
  obsTypes= st.auxTypes; 
  %obsData=st.auxData;
  obsData=st.auxData(:,iepoch,:);     % adjust to sorted epochs
else
  error('invalid item')
end

if ~isempty(opt.types)
   [~,idxTypes]=ismember(opt.types,obsTypes);
else
   idxTypes=1:numel(obsTypes);
end

% S-transform 

if ~strcmpi(opt.schranking,'none')
   % Select point reference mask
   pntRefMask = ismember(st.pntName,opt.refPoints);
   if exist(opt.unstablearea,'file')
      % Get polygon with unstable area
      unstable_area = load(opt.unstablearea);
      if isstruct(unstable_area)
         unstable_area = [ unstable_area.Lat(:) unstable_area.Lon(:) ];
      end
      pntRefMask = pntRefMask | ~inpolygon(st.pntCrd(:,2),st.pntCrd(:,1),unstable_area(:,2),unstable_area(:,1));
   end
   if strcmpi(opt.schranking,'min-velocities')
      fprintf('S-transform using min-velocities method: number of stable points: %d\n',sum(pntRefMask))
      [vel,t0,xoff,xref,omt]=stmvelocity(st, ...
         'refsystem',pntRefMask,'ignoreStochModel',true);
      xref=-1*xref(iepoch,:); % adjust sign and sorting
      obsData(:,:,idxTypes)=stmtransapply(obsData(:,:,idxTypes),xref(:,idxTypes),xoff(:,idxTypes));
   else
      % S-transform
      [xref,~,obsData(:,:,idxTypes)]=stmtrans(obsData(:,:,idxTypes),'method',opt.schranking,'pntRefMask',pntRefMask,'year',year);
   end
end

% Select stations for plotting

nstations=numel(pntName);
idxStation=1:nstations;
if ~isempty(opt.sort) 
   sortfield=strsplit(opt.sort,'(');
   if isfield(pntAttrib,sortfield{1})
      tmp=eval(sprintf('pntAttrib.%s',opt.sort));
      [~,idxStation]=sort(tmp);
   else
      warning('invalid pntAttrib in sort option, ignore.')
   end
end

% Prepare number of subplots

if numel(idxTypes) > opt.maxSubPlots
   error('Maximum number of items for plotting exceeded.')
end
if numel(idxTypes) > 1
   numsubplotspage = numel(idxTypes);
   maxstationsperpage= opt.maxTraces;
   numpages = min(ceil(nstations/maxstationsperpage), opt.maxPlots);
   numsubplots = numpages*numsubplotspage; 
   idxSubComponent = repmat(idxTypes(:),[numpages 1]);
   idxSubStation = [];
   for k=1:numpages
      idxSubStation = [idxSubStation ; repmat([ (k-1)*maxstationsperpage+1 min(nstations,k*maxstationsperpage) ],[ numsubplotspage 1]) ];
   end 
   numstationspersubplot=min(nstations,maxstationsperpage);
elseif numel(idxTypes) == 1
   numsubplots = ceil(nstations/opt.maxTraces); 
   numsubplotspage = min( numsubplots, opt.maxSubPlots); 
   numpages = min(ceil(numsubplots/numsubplotspage), opt.maxPlots);
   idxSubComponent = repmat(idxTypes(1),[numsubplots 1]);
   maxstationsperpage= opt.maxTraces * opt.maxSubPlots;
   numstationspersubplot=ceil(min(nstations,maxstationsperpage)/numsubplotspage);
   for k=1:numsubplots
      idxSubStation(k,1:2) = [ (k-1)*numstationspersubplot+1 min(k*numstationspersubplot,nstations) ];
   end
end

% numsubplots
% idxSubStation
% maxstationsperpage
% numstationspersubplot


% Prepare the plot, plot title and labels
sylabel= [];
for k=1:numel(obsTypes)
  sylabel{k}= [ obsTypes{k} ' [mm]'];
end

minyear=2999;
maxyear=0;

minyear=min(year);
maxyear=max(year);

% Make the plots
%
% - first loop over the three coordinate directions, creating a plot for each
% - then loop over the stations, and add a plot track for each station

if ~isempty(opt.saveplt)
   if ~exist(opt.saveplt,'dir')
       mkdir(opt.saveplt)
   end
   pdffile=fullfile(opt.saveplt,[ datasetId '_' item '.pdf']);
   fprintf('The plots will be saved to %s\n',pdffile);
   appendplt=false;
end

for kpage=1:numpages

    % Create the figure (in case of subplots)
    
    if opt.subplot
       hf=figure('Name',[datasetId '_' item '_' num2str(kpage)] , ...
              'NumberTitle','off',...
              'Position',[ opt.position(1) opt.position(2) min(numsubplotspage*opt.position(3),1860) opt.position(4) ]);
    end
    
    level=numstationspersubplot:-1:1;
    
    ksubplotpage=0;
    for ksubplot=(kpage-1)*numsubplotspage+1:min(numsubplots,kpage*numsubplotspage)
    
      icomponent=idxSubComponent(ksubplot);
    
      % Create a figure (in case of individual plots), or a subplot
      
      if opt.subplot
         ksubplotpage=ksubplotpage+1;
         subplot(1,numsubplotspage,ksubplotpage)
      else  
         hf(ksubplot)=figure('Name',[ datasetId ' ' obsTypes{icomponent} ], ...
                               'NumberTitle','off',...
                               'Position',[ opt.position(1)+(ksubplot-1)*(opt.position(3)+2) opt.position(2) opt.position(3) opt.position(4) ]);
      end
    
      % Add a track for each station
    
      ii=0;
      for i=idxSubStation(ksubplot,1):idxSubStation(ksubplot,2)
        ii=ii+1;
        istation=idxStation(i);
        data=obsData(istation,:,icomponent)*opt.scale;
        idxdata=find(~isnan(data));
        station=pntName{istation};
        offset=level(ii)*spacing(icomponent)-nanmean(data);
        if isempty(idxdata)
           text(mean(year),level(ii)*spacing(icomponent),[ station ' '],'HorizontalAlignment','right','Color','red') ;
           hold all;
           continue;
        end
        h=plot(year,data+offset,'LineWidth',1);
        hcol=get(h,'Color');
        hold all;
        line([year(idxdata(1)) ;year(idxdata(end))],[ level(ii)*spacing(icomponent) ; level(ii)*spacing(icomponent) ],'LineStyle',':','Color',hcol);
        if opt.marker
           plot(year(~isnan(data)),data(~isnan(data))+offset,'--o','Color',hcol);
        end
        % Plot symbol for events
        if isfield(pntAttrib,'events')
           events=strsplit(pntAttrib.events{istation},';')';
           numevents=numel(events);
        else
           numevents=0;
        end
        for k=1:numevents
           eventItems=split(events{k},',');
           eventType=eventItems{1};
           eventYear=date2dyear(datenum(eventItems{2},'yyyymmddTHHMMSS'));
           if ~ismember(eventType,opt.events), continue; end
           if numel(idxdata)<2, continue; end
           %ymark=level(i)*spacing(icomponent);
           ymark=interp1(year(idxdata),data(idxdata),eventYear,'linear')+offset;
           switch eventType
             case {'REC'}
                %plot(eventYear,ymark,'v','Color','y')
                text(eventYear,ymark,'R','Color',hcol,'FontSize',9,'HorizontalAlignment','center','VerticalAlignment','top')
             case {'ANT'}
                %plot(eventYear,ymark,'^','Color','r')
                text(eventYear,ymark,'A','Color',hcol,'FontSize',9,'HorizontalAlignment','center','VerticalAlignment','top')
             case {'DPL','STEP'}
                %plot(eventYear,ymark,'v','Color','k');
                text(eventYear,ymark,'^','Interpreter','none','Color',hcol,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','top')
             %otherwise
             %   plot(eventYear,ymark,'x','Color','g')
           end
           %text(eventYear,ymark,'|','Interpreter','none','Color',[.5 .5 .5],'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')
        end 
        text(year(idxdata(1)),level(ii)*spacing(icomponent),[ station ' '],'HorizontalAlignment','right','FontSize',opt.FontSize) ;
        minyear=min(minyear,floor(opt.xfmargin*year(idxdata(1)))/opt.xfmargin);
        maxyear=max(maxyear,ceil(opt.xfmargin*year(idxdata(end)))/opt.xfmargin);
      end
      
      % Do the title, plot labels, and beatify the plot
      
      title([ datasetId ' ' obsTypes{icomponent}  ],'Interpreter','none');
      ylabel(sylabel{icomponent});
      a=axis();
      a(1)=minyear-ceil(4*(maxyear-minyear)/15)/4;
      a(2)=maxyear;
      if rowmargin(1) < 0
         a(3)=min(a(3),(1-abs(rowmargin(1)))*spacing(icomponent));
      else
         a(3)=(1-rowmargin(1))*spacing(icomponent);
      end
      if rowmargin(2) < 0
         a(4)=max(a(4),(numstationspersubplot+abs(rowmargin(2)))*spacing(icomponent));
      else
         a(4)=(numstationspersubplot+rowmargin(2))*spacing(icomponent);
      end
      axis(a);
      set(gca,'Ylim',[a(3) a(4)],'Ytick',[1:numstationspersubplot]*spacing(icomponent));
      xlim=get(gca,'XLim');
      ylim=get(gca,'YLim');
      if xlim(2)-xlim(1) < 7
         xtick=xlim(1):1/12:xlim(2);
      else
         xtick=floor(xlim(1)):ceil(xlim(2));
      end
%       xticklabel=get(gca,'XTickLabel');
%       xticklabel
%       if iscellstr(xticklabel)
%          for k=1:size(xticklabel,1)
%            xticklabel{k}=sprintf('%6.1f',str2num(xticklabel{k}));
%          end
%       elseif ischar(xticklabel)
%          % fmt=sprintf('%%%d.%df',size(xticklabel,2),size(xticklabel,2)-5);
%          fmt='%6.1f';
%          for k=1:size(xticklabel,1)
%            xticklabelc(k,:)=sprintf(fmt,str2num(xticklabel(k,:)));
%          end
%          xticklabel=xticklabelc;
%       end
%       xticklabel
%      set(gca,'XTickLabel',xticklabel);
      ticklength=get(gca,'TickLength');
      set(gca,'YMinorTick','on');
      if opt.subplot
         %[ ((ksubplotpage-1)+opt.innerPosition(1))/numsubplotspage  opt.innerPosition(2) opt.innerPosition(3)/numsubplotspage  opt.innerPosition(4)]  
         set(gca,'InnerPosition', [ ((ksubplotpage-1)+opt.innerPosition(1))/numsubplotspage  opt.innerPosition(2) opt.innerPosition(3)/numsubplotspage  opt.innerPosition(4)] );
      else
         set(gca,'InnerPosition', opt.innerPosition);
      end
      cax1 = axes('Position', get(gca, 'Position'));
      cax2 = axes('Position', get(gca, 'Position'));
      set(cax1 ...
        ,'TickLength'       , ticklength/2 ...
        ,'YAxisLocation'    ,'left' ...
        ,'XAxisLocation'    , 'Top' ...
        ,'Box'              , 'off' ...
        ,'HitTest'          , 'off' ...
        ,'Box'              , 'off' ...
        ,'Color'            , 'none' ...
        ,'XTick'            , xtick ...
        ,'YTick'            , [] ...
        ,'XTickLabel'       , [] ...
        ,'XLim'             , xlim ...
        ,'YLim'             , ylim ...
        );
      set(cax2 ...
        ,'TickLength'       , ticklength/2 ...
        ,'YAxisLocation'    ,'left' ...
        ,'XAxisLocation'    , 'Bottom' ...
        ,'Box'              , 'off' ...
        ,'HitTest'          , 'off' ...
        ,'Box'              , 'off' ...
        ,'Color'            , 'none' ...
        ,'XTick'            , xtick ...
        ,'YTick'            , [] ...
        ,'XTickLabel'       , [] ...
        ,'XLim'             , xlim ...
        ,'YLim'             , ylim ...
        );

      if ~isempty(opt.saveplt) && ~opt.subplot 
          exportgraphics(hf(ksubplot),pdffile,'append',appendplt)
          appendplt=true;
      end

    end

    if ~isempty(opt.saveplt) && opt.subplot
       exportgraphics(hf,pdffile,'append',appendplt)
       appendplt=true;
    end

end

end

