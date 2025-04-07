function events=getEvents(eventfile,station,tfirst,tlast)
%getEvents   Get events from csv text file.
%   EVENTS=getEvent(EVENTFILE,STATION) reads the events for station STATION
%   from the EVENTFILE csv text file, and returns a structure array EVENTS
%   with the events for that station.
%
%   EVENTS=getEvent(EVENTFILE,STATION,TFIRST,TLAST) also trims the EVENTS
%   structure array to the range TFIRST - TLAST in decimal years.
%
%   The format for the csv EVENTFILE is
%
%     station,date,type,name,remarks
%     0645,2014-09-01T00:00Z,REC,LEICA GRX1200+GNSS  , 
%     0645,2015-06-24T00:00Z,REC,SEPT POLARX4        ,Inserted by HM 
%     0645,2014-09-01T00:00Z,ANT,LEIAR25.R4      LEIT, 
%     ...
%
%   See also prtEvents, trimEvents and addEvent. 

% Created:  19 September 2019 by Hans van der Marel
% Modified:

persistent loadedEventFile allEvents

if isempty(loadedEventFile) || ~strcmpi(eventfile,loadedEventFile)
   allEvents=readtable(eventfile,'format','%s %s %s %s %s');
   loadedEventFile=eventfile;
end

if nargin <=2 
  tfirst=-Inf;
  tlast=+Inf;
end

eventtable=table2struct(allEvents(ismember(allEvents.station,station),2:end));
for k=1:size(eventtable,1)
  events(k).year=date2dyear(eventtable(k).date,'yyyy-mm-ddTHH:MMZ');
  events(k).type=eventtable(k).type;
  events(k).name=eventtable(k).name;
end
events=trimEvents(events,tfirst,tlast,'station',station);

end

