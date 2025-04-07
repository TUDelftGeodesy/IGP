function prtEvents(ts,title,fieldname)
%prtEvents  Print the list of events for an event or time series structure.
%   prtEvents(EVENTS,TITLE) prints the EVENTS structure array with 
%   optional title string TITLE.
%
%   prtEvents(TS,TITLE) prints the EVENTS structure array within
%   the time series structure TS with title string TITLE.
%
%   prtEvents(TS,TITLE,FIELDNAME) prints the EVENTS structure array within
%   the time series structure TS with title string TITLE, but using
%   FIELDNAME as alternative fieldname in the structure.
%
%   See also getEvents, trimEvents and addEvent. 

% Created:  19 September 2019 by Hans van der Marel
% Modified:

if nargin < 3
   fieldname='events';
end

if nargin == 2, fprintf('Events %s:\n\n',title); end

if isfield(ts,fieldname)
  % is time series structure
  for k=1:numel(ts)
    for l=1:numel(ts(k).(fieldname))
      fprintf('%s   %10.3f   %5s   %s\n',ts(k).station,ts(k).(fieldname)(l).year,ts(k).(fieldname)(l).type,ts(k).(fieldname)(l).name)
    end
    fprintf('\n')
  end
else
  % must be event structure
  events=ts;
  for l=1:numel(events)
    fprintf('%10.3f   %5s   %s\n',events(l).year,events(l).type,events(l).name)
  end
 fprintf('\n') 
end

end

