function hf=tseriesplotcomp(ts,item,varargin)
%TSERIESPLOTCOMPONENT Plot timescomponents for multiple stations.
%   TSERIESPLOTCOMPONENT(TSERIES,ITEM) plots the timeseries North, East and
%   Up decomposed signal components, as specified in the cell array ITEM, 
%   for the stations in the structure array TSERIES.
%
%   TSERIESPLOTCOMPONENT(STATIONS,ITEM,OPTION,VALUE,...) accepts multiple
%   OPTION/VALUE pairs, valid option/value pairs are
%
%     Project       string       project name (in front of title)
%     Id            string       identifier of the fit (in brackets at end of title)
%     Spacing       scalar | 1x3 spacing between traces , default is 10 ([10 10 10])
%     RowMargin     scalar | 1x2 add rowmargin extra (blank) rows at top and bottom of plot (-1),
%                                scalar, or [bottom top], if negative, does auto scaling
%     Events        cell array   Events to plot, default ({'ANT' 'DPL' 'STEPS'}), i.e no 'REC' 
%     Subplot       logical      If true use subplots (default), if false create individual plots
%     Position      1x4 array    figure Position [ left bottom width height] in pixels ([50 50 560 870]) 
%                                for the leftmost plot (or subplot)
%     InnerPosition 1x4 array    axis InterPosition [ left bottom width height] in normalized units ([ 0.10 0.04 0.87 0.93]) 
%
%   HF=TSERIESPLOTCOMPONENT(...) returns the plot handles for the three
%   plots in the North, East and Up direction.
%
%   Examples:
%      tseriesplotcomponent(tseries,'raw');
%      tseriesplotcomponent(tseries,{'harmonic' 'tempi'},'id','it1');
%      tseriesplotcomponent(stations,{'harmonic' 'tempi'},'id','it1', ...
%                  'project','Groningen','spacing',10);
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

% Default options


opt.Project='';                % optional project name for title
opt.Id='';                     % optional id for title (will be overwritten by series id)
opt.Spacing=10;                % spacing between traces
opt.RowMargin=-1;              % add rowmargin extra (blank) rows at top and bottom of plot (-1)
opt.Events={'ANT' 'DPL' 'STEPS'};  % cell array with events to plot 
opt.Subplot=true;              % Use subplots (default), if false create individual plots
opt.Position=[50 50 560 870];   % figure Position [ left bottom width height] in pixels for the leftmost plot (or subplot)
%opt.InnerPosition=[ 0.13 0.11 0.775 0.815];   % axis InnerPostion [ left bottom width height] in normilized units 
opt.InnerPosition=[ 0.10 0.04 0.87 0.93];   % axis InnerPosition [ left bottom width height] in normalized units 

% Check arguments

if nargin < 2, error('too few arguments (specify at least stationames and item');end
if ~isstruct(ts)
    error('first argument must be a structure array of timeseries.') 
end
if nargin < 2
    item='raw';
end
item=cellstr(item);

% Check the option/value pairs

for k=1:2:numel(varargin)-1
   if isfield(opt,varargin{k})
      opt.(varargin{k})=varargin{k+1};
   else
      error(['Illegal option ' varargin{k} ])
   end
end

% Prepare the plot, plot title and labels

id=opt.Id;
project=opt.Project;
if isscalar(opt.Spacing)
  spacing=[opt.Spacing opt.Spacing opt.Spacing];
else
  spacing=opt.Spacing;
end
if isscalar(opt.RowMargin)
  rowmargin=[opt.RowMargin opt.RowMargin];
else
  rowmargin=opt.RowMargin;
end

itemlist=item{1};
for k=2:length(item)
  itemlist=[ itemlist ' + ' item{k} ];
end
if ~isempty(id)
  titlestr=[ ' (' itemlist ', id=' id ')'];
else
  titlestr=[ ' (' itemlist ')'];
end

for k=1:length(item)
  %if strcmp(item{k},'raw')
  if ismember(item,{'raw' 'nostep'})
    item{k}='neu';
  else
    item{k}= [ 'neu' item{k} ];
  end
end

remjump=false;
%if ismember(item,{'raw' 'fit'})
if ismember(item,{'nostep' 'fit'})
  remjump=true;
end

sylabel= { '\Delta N [mm]'  '\Delta E [mm]' '\Delta U [mm]' };
complabel= { 'North'  'East' 'Up' }; 

minyear=2999;
maxyear=0;

% Make the plots
%
% - first loop over the three coordinate directions, creating a plot for each
% - then loop over the stations, and add a plot track for each station

nstations=length(ts);          
level=nstations:-1:1;

% Create the figure (in case of subplots)

if opt.Subplot
   hf=figure('Name',[ project ' ' titlestr ], ...
          'NumberTitle','off',...
          'Position',[ opt.Position(1) opt.Position(2) 3*opt.Position(3) opt.Position(4) ]);
end

for icomponent=1:3

  % Create a figure (in case of individual plots), or a subplot
  
  if opt.Subplot
     subplot(1,3,icomponent)
  else  
     hf(icomponent)=figure('Name',[ project ' ' complabel{icomponent} ' ' titlestr ], ...
                           'NumberTitle','off',...
                           'Position',[ opt.Position(1)+(icomponent-1)*(opt.Position(3)+2) opt.Position(2) opt.Position(3) opt.Position(4) ]);
  end
  
  % Add a track for each station
  
  for i=1:nstations
    tseries=ts(i);
    station=tseries.station;
    if isfield(tseries,'year')
       year=tseries.year;
    else
       year=ymd2dyear(datevec(tseries.epoch));
    end
    data=tseries.(item{1});
    if remjump, data=data-tseries.neujumps; end
    for k=2:length(item)
       data=data+tseries.(item{k});
    end
    offset=level(i)*spacing(icomponent)-nanmean(data(:,icomponent))*1000;
    % insert NaN's at gaps
    idx=find(diff(year) > 1/360);
    n=length(idx);
    idx(n+1)=length(year);
    for k=n:-1:1
      year(idx(k)+k+1:idx(k+1)+k)=year(idx(k)+1:idx(k+1));
      data(idx(k)+k+1:idx(k+1)+k,:)=data(idx(k)+1:idx(k+1),:);
      year(idx(k)+k)=nan;
      data(idx(k)+k,:)=nan;
    end
    h=plot(year,data(:,icomponent)*1000+offset,'LineWidth',1);
    hcol=get(h,'Color');
    line([year(1) ;year(end)],[ level(i)*spacing(icomponent) ; level(i)*spacing(icomponent) ],'LineStyle',':','Color',hcol);
    hold all;
    % Plot symbol for events
    events=tseries.events;
    for k=1:numel(events)
       if ~ismember(events(k).type,opt.Events), continue; end
       %ymark=level(i)*spacing(icomponent);
       ymark=interp1(year(~isnan(year)),data(~isnan(year),icomponent)*1000,events(k).year,'linear')+offset;
       switch events(k).type
         case {'REC'}
            %plot(events(k).year,ymark,'v','Color','y')
            text(events(k).year,ymark,'R','Color',hcol,'FontSize',9,'HorizontalAlignment','center','VerticalAlignment','top')
         case {'ANT'}
            %plot(events(k).year,ymark,'^','Color','r')
            text(events(k).year,ymark,'A','Color',hcol,'FontSize',9,'HorizontalAlignment','center','VerticalAlignment','top')
         case {'DPL','STEP'}
            %plot(events(k).year,ymark,'v','Color','k');
            text(events(k).year,ymark,'^','Interpreter','none','Color',hcol,'FontSize',14,'HorizontalAlignment','center','VerticalAlignment','top')
         %otherwise
         %   plot(events(k).year,ymark,'x','Color','g')
       end
       %text(events(k).year,ymark,'|','Interpreter','none','Color',[.5 .5 .5],'FontSize',12,'HorizontalAlignment','center','VerticalAlignment','middle')
    end 
    text(year(1),level(i)*spacing(icomponent),[ station ' '],'HorizontalAlignment','right') ;
    minyear=min(minyear,floor(4*min(year))/4);
    maxyear=max(maxyear,ceil(4*max(year))/4);
  end
  
  % Do the title, plot labels, and beatify the plot
  
  title([ project ' ' complabel{icomponent} ' ' titlestr],'Interpreter','none');
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
     a(4)=max(a(4),(nstations+abs(rowmargin(2)))*spacing(icomponent));
  else
     a(4)=(nstations+rowmargin(2))*spacing(icomponent);
  end
  axis(a);
  set(gca,'Ylim',[a(3) a(4)],'Ytick',[1:nstations]*spacing(icomponent));
  xlim=get(gca,'XLim');
  ylim=get(gca,'YLim');
  xtick=xlim(1):1/12:xlim(2);
  ticklength=get(gca,'TickLength');
  xticklabel=get(gca,'XTickLabel');
  if iscellstr(xticklabel)
     for k=1:size(xticklabel,1)
       xticklabel{k}=sprintf('%6.1f',str2num(xticklabel{k}));
     end
  elseif ischar(xticklabel)
     % fmt=sprintf('%%%d.%df',size(xticklabel,2),size(xticklabel,2)-5);
     fmt='%6.1f';
     for k=1:size(xticklabel,1)
       xticklabelc(k,:)=sprintf(fmt,str2num(xticklabel(k,:)));
     end
     xticklabel=xticklabelc;
  end
  set(gca,'XTickLabel',xticklabel);
  set(gca,'YMinorTick','on');
  if opt.Subplot
     set(gca,'InnerPosition', [ (icomponent-1)/3+opt.InnerPosition(1)/3  opt.InnerPosition(2) opt.InnerPosition(3)/3  opt.InnerPosition(4)] );
  else
     set(gca,'InnerPosition', opt.InnerPosition);
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
  
end

end