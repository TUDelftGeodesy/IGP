function pltseries1(year,dneu,sneu,rcvrtypes,ircvr,anttypes,iant,limits,sigopt)
%PLTSERIES1   Elementary timeseries plot with annotations

% sigopt=2 - plot two sigma fill (+- 2*sneu) around the plotdata (dneu)
% sigopt=1 - plot one sigma fill (+- sneu) around zero
% sigopt=0 - plot no fill at all

dd=3.5/365;
linetype='b-';

if isempty(dneu)
  dneu=zeros(size(year));
  linetype='c-';
end

% Get axis dimension

plot(year,dneu*1000,linetype)
a=axis;
if nargin > 7
   a=[ a(1:2) limits];
   axis(a);
end
if nargin < 9
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
  
% Report receiver and antenna changes

for k=1:length(ircvr)
   ip=ircvr(k);
   line([ ip-dd ; ip-dd], a(3:4),'LineStyle',':','Color','y')
   if ~isempty(rcvrtypes)
      text(ip,a(3)+(a(4)-a(3))*.8,rcvrtypes(k,:),'Interpreter','none','FontSize',7)
   end
end
for k=1:length(iant)
   ip=iant(k);
   line([ ip-dd ; ip-dd], a(3:4),'LineStyle',':','Color','r')
   if ~isempty(anttypes)
      text(ip,a(3)+(a(4)-a(3))*.9,anttypes(k,:),'Interpreter','none','FontSize',7) %'Color','r')
   end
end

axis(a)

end
