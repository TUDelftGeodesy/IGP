function prtSteps(tsfit,titlestr)
%prtSteps   Print the estimated steps for a time series structure array.
%   prtSteps(TSFIT) prints the estimated steps for the time series
%   structure array TSFIT.
%
%   prtSteps(TSFIT,TITLESTR) uses TITLESTR instead of time series id as
%   title.

if nargin < 1
   error('Function expects a time series structure as input.')
end
if ~isfield(tsfit,'parindex')
   error('Time series structure does not contain estimated parameters.')
end
if nargin < 2
   titlestr=char(unique({tsfit.id}));
end

fprintf('Estimated steps %s:\n\n',titlestr)
fprintf(' sid     dyear event      N[mm]     E[mm]     U[mm]  sN[mm] sE[mm] sU[mm]\n');
fprintf('----  --------  ----  --------- --------- ---------  ------ ------ ------\n');
for k=1:numel(tsfit)
  station=tsfit(k).station;
  eventtypes={ tsfit(k).events(:).type };
  eventfilter=ismember(eventtypes,{'ANT','DPL','STEP','REF'});
  eventtypes=eventtypes(eventfilter);
  eventyears=[ tsfit(k).events(:).year ];
  ll=0;
  for l=tsfit(k).parindex(3)+1:tsfit(k).parindex(4)
    ll=ll+1;
    tjump=tsfit(k).tjump(ll);  
    idxevent=find(abs(eventyears-tjump) < 1/365);
    step=[ tsfit(k).xlat(l) tsfit(k).xlon(l) tsfit(k).xrad(l)  ];
    sigstep=sqrt([ tsfit(k).qxlat(l,l) tsfit(k).qxlon(l,l) tsfit(k).qxrad(l,l)].*tsfit(k).omt); 
    fprintf('%s  %8.3f  %4s  %9.2f %9.2f %9.2f  %6.2f %6.2f %6.2f\n',station,tjump,eventtypes{idxevent},step*1000,sigstep*1000);
  end
end

end