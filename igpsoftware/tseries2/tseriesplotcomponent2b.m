function hf=tseriesplotcomponent2b(ts1,ts2,item,project,spacing)
%TSERIESPLOTCOMPONENT2 Plot dual timescomponents for multiple stations.
%   TSERIESPLOTCOMPONENT(TSERIES1,TSERIES2,ITEM) plots the timeseries components
%   in the cell array ITEM for the stations in the structure array TSERIES1
%   and TSERIES2 over each other. TSERIES1 is plotter over TSERIES2.
%
%   TSERIESPLOTCOMPONENT(...,PROJECT,SPACING) provides additional project
%   name PROJECT and spacing SPACING between the individual lines. Default
%   for spacing is 10 ([10 10 10]).
%
%   Examples:
%      tseriesplotcomponent(tseries1,tseries2,'raw');
%      tseriesplotcomponent(tseries1,tseries2,{'harmonic' 'tempi'});
%      tseriesplotcomponent(tseries1,tseries2,{'harmonic' 'tempi'},'Ameland',10);
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2019.

id='';

% Check arguments

if nargin < 2, error('too few arguments (specify at least stationames and item');end
if ~isstruct(ts1) && ~isstruct(ts2)
    error('first and second argument must be structure arrays of timeseries.') 
end
if nargin < 3
    item='raw';
end
item=cellstr(item);
if nargin < 4
   project='';
end
if nargin < 5
  spacing=10;
end
if isscalar(spacing)
  spacing=[spacing spacing spacing];
end

saveplots=false;

% Get station names, set labels and items

stations=unique([{ ts1.station } { ts2.station }],'stable'); 

nstations=length(stations);          
level=1:nstations;

sylabel= { '\Delta N [mm]'  '\Delta E [mm]' '\Delta U [mm]' }; 
complabel= { 'North'  'East' 'Up' }; 

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

minyear=2999;
maxyear=0;

remjump=false;
%if ismember(item,{'raw' 'fit'})
if ismember(item,{'nostep' 'fit'})
  remjump=true;
end

% Get colors

col=get(groot,'defaultAxesColorOrder');
nc=size(col,1);
nc2=floor(nc/2);

% Make a plot of each component

for icomponent=1:3
  hf(icomponent)=figure;
  for i=1:nstations
    station=stations{i};
    y1=+Inf;
    y2=-Inf;
    year1=[];
    data1=[];
    if any(ismember({ts2.station},station))
      tseries=ts2(ismember({ts2.station},station));
      if isfield(tseries,'year'), year=tseries.year; else year=ymd2dyear(datevec(tseries.epoch)); end
      data=tseries.(item{1});
      if remjump, data=data-tseries.neujumps; end
      for k=2:length(item)
        data=data+tseries.(item{k});
      end
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
      offset=level(i)*spacing(icomponent)-nanmean(data(end-10:end,icomponent))*1000;
      hcol=col(mod(i,7)+1,:);
      h=plot(year,data(:,icomponent)*1000+offset,'LineWidth',1,'Color',brighten(hcol,.5));
      %hcol=get(h,'Color');
      hold all; 
      y1=min(y1,year(1));
      y2=max(y2,year(end));
      year1=year;
      data1=data;
    end      
    if any(ismember({ts1.station},station))
      tseries=ts1(ismember({ts1.station},station));
      if isfield(tseries,'year'), year=tseries.year; else year=ymd2dyear(datevec(tseries.epoch)); end
      data=tseries.(item{1});
      if remjump, data=data-tseries.neujumps; end
      for k=2:length(item)
        data=data+tseries.(item{k});
      end
      % insert NaN's at gaps
      idx=find(diff(year) > 1.5/365);
      n=length(idx);
      idx(n+1)=length(year);
      for k=n:-1:1
        year(idx(k)+k+1:idx(k+1)+k)=year(idx(k)+1:idx(k+1));
        data(idx(k)+k+1:idx(k+1)+k,:)=data(idx(k)+1:idx(k+1),:);
        year(idx(k)+k)=nan;
        data(idx(k)+k,:)=nan;
      end
      if ~isempty(year1)
         [~,ia,ib]=intersect(round(dyear2date(year1)),round(dyear2date(year)));
         offset=level(i)*spacing(icomponent)-nanmean(data1(end-10:end,icomponent))*1000+nanmean(data1(ia,icomponent)-data(ib,icomponent))*1000;         
      else
         offset=level(i)*spacing(icomponent)-nanmean(data(end-10:end,icomponent))*1000;
      end
      hcol=col(mod(i,7)+1,:);
      h=plot(year,data(:,icomponent)*1000+offset,'LineWidth',1,'Color',brighten(hcol,-.5));
      %hcol=get(h,'Color');
      hold all;
      y1=min(y1,year(1));
      y2=max(y2,year(end));
    end      
    line([y1 ;y2],[ level(i)*spacing(icomponent) ; level(i)*spacing(icomponent) ],'LineStyle',':','Color',hcol);
%   text(y1,level(i)*spacing(icomponent),[ station ' '],'HorizontalAlignment','right') ;
%   text(y1,level(i)*spacing(icomponent)+nanmean(data(1:10,icomponent))*1000-nanmean(data(end-10:end,icomponent))*1000,[ station ' '],'HorizontalAlignment','right') ;
    minyear=min(minyear,floor(4*y1)/4);
    maxyear=max(maxyear,ceil(4*y2)/4);
  end
  title([ project ' ' complabel{icomponent} ' ' titlestr],'Interpreter','none');
  ylabel(sylabel{icomponent});
  a=axis();
  %a(1)=minyear-ceil(4*(maxyear-minyear)/15)/4;
  a(1)=minyear;
  a(2)=maxyear;
  a(3)=max(min(a(3),0),-spacing(icomponent));
  a(4)=(nstations+1.5)*spacing(icomponent);
  axis(a);
  set(gca,'Ylim',[a(3) a(4)],'Ytick',[0:10:a(4)]);
  for i=1:nstations
    station=stations{i};
    text(maxyear,level(i)*spacing(icomponent),[' ' station],'HorizontalAlignment','left') ;
  end
  
  xlim=get(gca,'XLim');
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
    );
  
end

if saveplots
  plotdir='plots';
  filename=['all_' strrep(itemlist,' ','') id ]; 
  print(hf(1),'-dpng','-r300',fullfile(plotdir,[ filename '_lat.png']));
  print(hf(2),'-dpng','-r300',fullfile(plotdir,[ filename '_lon.png']));
  print(hf(3),'-dpng','-r300',fullfile(plotdir,[ filename '_hgt.png']));
end

end