function pltspectrum(t,y,ptitle,doplomb)
%PLTSPECTRUM  Plot GPS Lomb-Scargle Periodogram

% Compute Lomb-Scargle periodogram using Matlab function plomb
%
% Older versions of Matlab don't support plomb, therefore the Matlab 2014a
% version plomb.m has been copied to the current directory, with lines 80-85
% commented out, and the message catalog plomb.xml is copied to 
% resources/signal/en/plomb.xml. 

if nargin < 4 
  doplomb=true;
end
if doplomb
  [s,f]=plomb([y(:,1)*10 y(:,2)*100 y(:,3)*1000],t);
else
  f=t;
  s=y;
end

ylow=1e-8;
yup=100;

figure
plot(f,s)
set(gca,'XScale','log')
set(gca,'YScale','log')
xlabel('Frequency (cpy)')
ylabel('Power/Frequency (mm^2/cpy)')
set(gca,'Xlim',[0.01 200 ])
set(gca,'YLim',[ ylow yup])
legend('N','E','U')
title(ptitle,'Interpreter','none')

ylow=ylow*2;
yup=yup/4;
line( [1 1], [ylow yup],'LineStyle',':','Color','k')
line( [2 2], [ylow yup],'LineStyle',':','Color','k')
text(1,yup*2,'1','HorizontalAlignment','center')
text(2,yup*2,'2','HorizontalAlignment','center')

line( [1.04 1.04], [ylow yup],'LineStyle',':','Color','k')
line( [2.08 2.08], [ylow yup],'LineStyle',':','Color','k')
line( [3.12 3.12], [ylow yup],'LineStyle',':','Color','k')
line( [4.18 4.18], [ylow yup],'LineStyle',':','Color','k')
line( [5.20 5.20], [ylow yup],'LineStyle',':','Color','k')
line( [6.24 6.24], [ylow yup],'LineStyle',':','Color','k')
text(3.12,yup*2,'3.12','HorizontalAlignment','center')
text(4.16,yup*2,'4.16','HorizontalAlignment','center')
text(5.20,yup*2,'5.20','HorizontalAlignment','center')
text(6.24,yup*2,'6.24','HorizontalAlignment','center')

line( [24.76 24.76], [ylow yup],'LineStyle',':','Color','k')
text(24.76,yup*2,'14.75d','HorizontalAlignment','center')

line(f,1e-1*f.^-1,'LineStyle',':','Color','k')
line(f,1e-3*f.^-1,'LineStyle',':','Color','k')
line(f,1e-5*f.^-1,'LineStyle',':','Color','k')
