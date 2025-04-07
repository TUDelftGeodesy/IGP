function [events,idx]=addEvent(ts,station,year,type,name)
%addEvent   Add event to list of events from a time series structure.
%   [EVENTS,IDX]=addEvent(TS,STATION,YEAR,TYPE,NAME) adds the event at
%   decimal year YEAR, type TYPE and name NAME to the events structure
%   in the time series structure array TS for station with id STATION.
%   Output is the modified events structure EVENTS for TS(IDX). The
%   original time series structure is not modified.
%
%   See also getEvents, prtEvents and trimEvents. 

% Created:  19 September 2019 by Hans van der Marel
% Modified:

[stationname,idx]=tsSelect(ts,station);
if numel(idx) ~= 1, disp(['station ' station ' not found']); events=[]; idx=[]; return; end

events=ts(idx).events;

[~,k]=sort([events.year]);
events=events(k);

k1=find([events.year] <= year);
knew=k1(end)+1;
events(knew+1:end+1)=events(knew:end);
events(knew).year=year;
events(knew).type=type;
events(knew).name=name;


fprintf('New event added to %s:\n\n',char(stationname));
for l=1:numel(events)
  if l==knew, flag='** NEW EVENT **'; else flag=''; end
  fprintf('%s   %10.3f   %5s   %20s  %s\n',char(stationname),events(l).year,events(l).type,events(l).name,flag)
end
fprintf('\n')

end
