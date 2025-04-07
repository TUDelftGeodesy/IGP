function pltseries1b(year,dneu,sneu,events,limits,sigopt)
%PLTSERIES1B  Elementary timeseries plot with annotations
%  PLTSERIES1B plots one component of the time series with a cyan fill
%  around as indication of the error. YEAR is an array with the decimal 
%  year, DNEU an array with the timeseries and SNEU an array representing
%  the error (one sigma). EVENTS is a structure array with equipment 
%  information, and other events, which will be used for annotations.  
%  LIMITS is an array with [ lower upper ] minimum limits, and SIGOPT 
%  specifies options for plotting the errors. 
%
%  All inputs are in [m] (where appropriate), but are converted into [mm]
%  in the plot.
%
%  EVENTS is a structure array with at least the following fields
%    events(k).year   start of the k'th event in decimal years
%    events(k).type   type of event: REC=receiver change, ANT=antenna
%                     change, or others.
%    events(k).name   name of the receiver or antenna (only for REC or ANT)
%
%  SIGOPT provides the following options 
%    sigopt=2    plot two sigma fill (+- 2*sneu) around the plotdata (dneu)
%    sigopt=1    plot one sigma fill (+- sneu) around zero
%    sigopt=0    plot no fill at all
%
%  If NEU is empty only the cyan fill with errors is plotted. If SNEU is
%  empty no fill is plotted.
%
%  This function is a lower level plot function that is usually only
%  called by other functions.

dd=3.5/365;
linetype='b-';

if isempty(dneu)
  dneu=zeros(size(year));
  linetype='c-';
end

% Get axis dimension

plot(year,dneu*1000,linetype)
a=axis;
if nargin > 4
   a=[ a(1:2) limits];
   axis(a);
else
   limits=[-Inf Inf];
end
if nargin < 6
   sigopt=2;
end

% Plot the data

dtmean=median(diff(year))*1.9;
iseg=[ 0 ; find(diff(year) > dtmean) ; length(year) ]; 
for i=1:length(iseg)-1;
  i1=iseg(i)+1;
  i2=iseg(i+1);  
  if ~isempty(sneu)
    if sigopt == 2
      fill([year(i1:i2);year(i2:-1:i1)],[dneu(i1:i2)*1000+sneu(i1:i2)*2000 ; dneu(i2:-1:i1)*1000-sneu(i2:-1:i1)*2000],'c','EdgeColor','c');
      hold on
    elseif sigopt == 1
      fill([year(i1:i2);year(i2:-1:i1)],[sneu(i1:i2)*1000 ; -sneu(i2:-1:i1)*1000],'c','EdgeColor','c');
      hold on
    end
  end
  plot(year(i1:i2),dneu(i1:i2)*1000,linetype)    
  hold on
end

if min(dneu)*1000/limits(1) > 1.1 || max(dneu)*1000/limits(2) > 1
   a=[ a(1:2) round(min(dneu)*1000) round(max(dneu)*1000) ] ;  
end
  
% Report receiver and antenna changes, and other events

for k=1:numel(events)
   ip=events(k).year;
   switch events(k).type
     case {'REC'}
       line([ ip-dd ; ip-dd], a(3:4),'LineStyle',':','Color','y')
       text(ip,a(3)+(a(4)-a(3))*.8,events(k).name,'Interpreter','none','FontSize',7)
     case {'ANT'}
       line([ ip-dd ; ip-dd], a(3:4),'LineStyle','--','Color','r')
       text(ip,a(3)+(a(4)-a(3))*.9,events(k).name,'Interpreter','none','FontSize',7) %'Color','r')
     case {'DPL'}
       line([ ip-dd ; ip-dd], a(3:4),'LineStyle','--','Color','r')
       %text(ip,a(3)+(a(4)-a(3))*.9,'DPL','Interpreter','none','FontSize',7) %'Color','r')
     case {'STEP'}
       line([ ip-dd ; ip-dd], a(3:4),'LineStyle','--','Color','r')
       text(ip,a(3)+(a(4)-a(3))*.9,'STEP','Interpreter','none','FontSize',7) %'Color','r')
     otherwise
       line([ ip-dd ; ip-dd], a(3:4),'LineStyle',':','Color','m')
   end
end

axis(a)

end
