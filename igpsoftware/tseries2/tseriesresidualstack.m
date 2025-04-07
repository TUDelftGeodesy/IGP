function [rstack,hfig]=tseriesresidualstack(stations,id,project,doplots,saveplotdir)
 %TSERIESRESIDUALSTACK  Compute and plot residual stack from fitted timeseries.
%   TSERIESRESIDUALSTACK(TSERIES,ID) computes and plots a residual stack
%   from the timeseries in the structure array TSERIES. ID is an optional 
%   for the output stack identifier which is saved as Mat file
%   rstack_id.mat. If ID is missing the residual stack is not saved to
%   file.
%
%   TSERIESRESIDUALSTACK(STATIONS,ID) computes and plots a residual stack
%   from the timeseries for the stations in the cell array STATIONS. 
%   ID is the optional identifier of the fit (default '_fit'). The residual
%   stack is saved to the mat file rstack_ID.mat. The timeseries
%   data must be stored in an mat-file STATION_ID.mat (_ is part of the
%   ID).
%
%   TSERIESRESIDUALSTACK(...,PROJECT) provides additional project
%   name PROJECT.
%
%   [RSTACK,HFIG]=TSERIESRESIDUALSTACK(...) also returns the residual
%   stack and figure handles.
%
%   Examples:
%      tseriesresidualstack(tseries,'_fit');
%      tseriesresidualstack(stations,'_fit');
%
%   (c) Hans van der Marel, Delft University of Technology, 2015-2016, 2021.

% Check arguments

if nargin < 1, error('too few arguments (specify at least stationames or timeseries.');end
if ~iscell(stations) && ~isstruct(stations)
    error('first argument must be a cell array with station names or a structure array of timeseries.') 
end
if nargin < 2
    item='_fit';
end
if nargin < 3
   project='';
end
if nargin < 4
   doplots=0;
end
if nargin < 5
   saveplotdir=0;
end

if ~isempty(project)
  titlestr=[ project  '(Stacked Residuals, id=' id ')' ];
  plotname=[ project '_rstack' ];
else
  titlestr=[ 'Stacked Residuals (id=' id ')' ];
  plotname='rstack'; 
end

saveplots=false;

% Compute the stack

nstations=length(stations); 
if iscell(stations)
  stationlist=stations;
elseif isstruct(stations)
  stationlist={ stations.station };
end

year=[];
data=[];
ghour=[];
sindex=zeros(nstations+1,1);
for i=1:nstations
  if iscell(stations)
     station=stations{i};
     load([ station id '.mat']);
  elseif isstruct(stations)
     tseries=stations(i);
     station=tseries.station;
  end
  year=[year ; tseries.year];
  ghour=[ ghour ; round((datenum(dyear2ymd(tseries.year))-datenum([2013 1 1]))*24)];
  data=[data ; tseries.neuresidual];
  sindex(i+1)=length(year);
end

[rmean,rstd,rcount,shour]=grpstats(data,ghour,{'mean', 'std', 'numel', 'gname'});
repoch=datenum([2013 1 1])+str2num(char(shour))./24;
ryear=ymd2dyear(datevec(repoch));

gday=floor(ghour/24);
[rmeanday,rstdday,rcountday,sday]=grpstats(data,gday,{'mean', 'std', 'numel', 'gname'});
repochday=datenum([2013 1 1])+str2num(char(sday));
ryearday=ymd2dyear(datevec(repochday));

rstack.titlestr=titlestr;
rstack.stations=stationlist;
rstack.sindex=sindex;
rstack.interval=median(diff(ghour));

rstack.year=year;
rstack.data=data;
rstack.ghour=ghour;

rstack.ryear=ryear;
rstack.rmean=rmean;
rstack.rstd=rstd;
rstack.rcount=rcount;
rstack.rhour=str2num(char(shour));

rstack.ryearday=ryearday;
rstack.rmeanday=rmeanday;
rstack.rstdday=rstdday;
rstack.rcountday=rcountday;
rstack.rday=str2num(char(sday));

% Save the residual stack to a mat file

%if nargin > 1
%  save(['rstack' id '.mat'],'rstack')
%end

% Plot the residual stack

if doplots > 0

    hfig(1)=figure('Name',[ plotname '_daily_' id ],'NumberTitle','off');
    subplot(3,1,1)
    pltseries1(ryearday,rmeanday(:,1),rstdday(:,1),[],[],[],[],[-4 4],1)
    h = findobj(gca,'Type','patch');
    h1 = findobj(gca,'Type','line','-and','Color','blue');
    a=axis;
    axis([ floor(5*min(ryearday))/5 ceil(5*max(ryearday))/5 a(3:4)]);
    ylabel('\Delta N [mm]')
    title(titlestr,'Interpreter','none')
    hl=legend([ h1(1) h(1)],'Mean (daily)','Std (daily)','Location','SouthWest','Orientation','horizontal');
    plabel=get(hl,'Position');
    plabel(1)=0.5-plabel(3)/2;
    plabel(2)=0.02;
    set(hl,'Position',plabel);
    subplot(3,1,2)
    pltseries1(ryearday,rmeanday(:,2),rstdday(:,2),[],[],[],[],[-4 4],1)
    a=axis;
    axis([ floor(5*min(ryearday))/5 ceil(5*max(ryearday))/5 a(3:4)]);
    ylabel('\Delta E [mm]')
    subplot(3,1,3)
    pltseries1(ryearday,rmeanday(:,3),rstdday(:,3),[],[],[],[],[-4 4],1)
    a=axis;
    axis([ floor(5*min(ryearday))/5 ceil(5*max(ryearday))/5 a(3:4)]);
    ylabel('\Delta U [mm]')

    ymd=dyear2ymd(year);
    ghourinmonth=round((datenum(ymd)-datenum([ymd(:,1) ymd(:,2) ones(size(ymd(:,3)))]))*24);
    [rmeanmonth,rstdmonth,smonth]=grpstats(data,ghourinmonth,{'mean', 'std', 'gname'});
    rdaymonth=str2num(char(smonth))./24;

    hfig(2)=figure('Name',[ plotname '_dayinmonth_' id ],'NumberTitle','off');
    subplot(3,1,1)
    pltseries1(rdaymonth,rmeanmonth(:,1),rstdmonth(:,1),[],[],[],[],[-2 2],1)
    set(gca,'XLim',[-0.5 32])
    h = findobj(gca,'Type','patch');
    h1 = findobj(gca,'Type','line','-and','Color','blue');
    ylabel('\Delta N [mm]')
    title([ titlestr ' - Monthly'],'Interpreter','none')
    hl=legend([ h1(1) h(1)],'Mean','Std','Location','SouthWest','Orientation','horizontal');
    plabel=get(hl,'Position');
    plabel(1)=0.5-plabel(3)/2;
    plabel(2)=0.02;
    set(hl,'Position',plabel);
    subplot(3,1,2)
    pltseries1(rdaymonth,rmeanmonth(:,2),rstdmonth(:,2),[],[],[],[],[-2 2],1)
    set(gca,'XLim',[-0.5 32])
    ylabel('\Delta E [mm]')
    subplot(3,1,3)
    pltseries1(rdaymonth,rmeanmonth(:,3),rstdmonth(:,3),[],[],[],[],[-2 2],1)
    set(gca,'XLim',[-0.5 32])
    ylabel('\Delta U [mm]')

    hfig(3)=figure('Name',[ plotname '_numel_' id ],'NumberTitle','off');
    subplot(2,1,1)
    plot(rstack.ryearday,rstack.rcountday)
    ylabel('numel')
    title([ titlestr ' - Number of elements'],'Interpreter','none')
    subplot(2,1,2)
    hist(rstack.rcountday(:,1),1:max(rstack.rcountday(:,1)));
    ylabel('count')
    xlabel('numel')

    if rstack.interval < 12

        hfig(4)=figure('Name',[ plotname '_hourly_' id ],'NumberTitle','off');
        subplot(3,1,1)
        pltseries1(ryear,rmean(:,1),rstd(:,1),[],[],[],[],[-4 4],1)
        h = findobj(gca,'Type','patch');
        h1 = findobj(gca,'Type','line','-and','Color','blue');
        h2=plot(ryearday,rmeanday(:,1)*1000,'g-');
        a=axis;
        axis([ floor(5*min(ryear))/5 ceil(5*max(ryear))/5 a(3:4)]);
        ylabel('\Delta N [mm]')
        title(titlestr,'Interpreter','none')
        hl=legend([ h1(1) h2(1) h(1)],'Mean (hourly)','Mean (daily)','Std (hourly)','Location','SouthWest','Orientation','horizontal');
        plabel=get(hl,'Position');
        plabel(1)=0.5-plabel(3)/2;
        plabel(2)=0.02;
        set(hl,'Position',plabel);
        subplot(3,1,2)
        pltseries1(ryear,rmean(:,2),rstd(:,2),[],[],[],[],[-4 4],1)
        plot(ryearday,rmeanday(:,2)*1000,'g-');
        a=axis;
        axis([ floor(5*min(ryear))/5 ceil(5*max(ryear))/5 a(3:4)]);
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        pltseries1(ryear,rmean(:,3),rstd(:,3),[],[],[],[],[-4 4],1)
        plot(ryearday,rmeanday(:,3)*1000,'g-');
        a=axis;
        axis([ floor(5*min(ryear))/5 ceil(5*max(ryear))/5 a(3:4)]);
        ylabel('\Delta U [mm]')

        ghourinday=round((datenum(ymd)-datenum([ymd(:,1) ymd(:,2) ymd(:,3)]))*24);
        [rmeandiurnal,rstddiurnal,sdiurnal]=grpstats(data,ghourinday,{'mean', 'std', 'gname'});
        rdiurnal=str2num(char(sdiurnal));

        hfig(5)=figure('Name',[ plotname '_diurnal1_' id ],'NumberTitle','off');
        subplot(3,1,1)
        pltseries1(rdiurnal,rmeandiurnal(:,1),rstddiurnal(:,1),[],[],[],[],[-2 2],1)
        h = findobj(gca,'Type','patch');
        h1 = findobj(gca,'Type','line','-and','Color','blue');
        ylabel('\Delta N [mm]')
        title([ titlestr ' - Diurnal'],'Interpreter','none')
        hl=legend([ h1(1) h(1)],'Mean','Std','Location','SouthWest','Orientation','horizontal');
        plabel=get(hl,'Position');
        plabel(1)=0.5-plabel(3)/2;
        plabel(2)=0.02;
        set(hl,'Position',plabel);
        subplot(3,1,2)
        pltseries1(rdiurnal,rmeandiurnal(:,2),rstddiurnal(:,2),[],[],[],[],[-2 2],1)
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        pltseries1(rdiurnal,rmeandiurnal(:,3),rstddiurnal(:,3),[],[],[],[],[-2 2],1)
        ylabel('\Delta U [mm]')

        hfig(6)=figure('Name',[ plotname '_diurnal2_' id ],'NumberTitle','off');
        subplot(3,1,1)
        boxplot(data(:,1)*1000,ghourinday)
        set(gca,'YLim',[-5 5])
        ylabel('\Delta E [mm]')
        title([ titlestr ' - Diurnal'],'Interpreter','none')
        subplot(3,1,2)
        boxplot(data(:,2)*1000,ghourinday)
        set(gca,'YLim',[-5 5])
        ylabel('\Delta E [mm]')
        subplot(3,1,3)
        boxplot(data(:,3)*1000,ghourinday)
        set(gca,'YLim',[-5 5])
        ylabel('\Delta U [mm]')

    end

    % Save the plots

    if ~isempty(saveplotdir)
      print(hfig(1),fullfile(saveplotdir,[ plotname '-daily' id '.png']),'-dpng','-r600')
      print(hfig(2),fullfile(saveplotdir,[ plotname '-dayinmonth' id '.png']),'-dpng','-r600')
      if rstack.interval > 12 
        print(hfig(3),fullfile(saveplotdir,[ plotname '-hourly' id '.png']),'-dpng','-r600')
        print(hfig(4),fullfile(saveplotdir,[ plotname '-diurnal1' id '.png']),'-dpng','-r600')
        print(hfig(5),fullfile(saveplotdir,[ plotname '-diurnal2' id '.png']),'-dpng','-r600')
      end
    end
end

end