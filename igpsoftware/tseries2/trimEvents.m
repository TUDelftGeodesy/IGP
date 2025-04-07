function updevents=trimEvents(events,tfirst,tlast,varargin)
%trimEvents   Trim the list of events for a given data period.
%   EVENTS=trimEvents(EVENTS,TFIRST,TLAST) trims the EVENTS structure
%   array to the data period between decimal year TFIRST and TLAST.
%
%   See also getEvents, prtEvents and addEvent. 

% Created:  19 September 2019 by Hans van der Marel
% Modified: 28 September 2019 by Hans van der Marel
%              - bugfix regarding first ANT and REC to tfirst
%              - if tfirst is before first ANT and REC entry, reset these
%                entries to tfirst (assuming it is the same antenna...)

% Check input arguments

if nargin <=1
   tfirst=-Inf;
   tlast=+Inf;
end

% Check options

opt.verbose=0;
opt.station='';
for k=1:2:numel(varargin)
   opt.(varargin{k})=varargin{k+1};    
end
if ~strcmpi(opt.station,''), opt.verbose=1; end

% Sort the events on type and decimal year (in that order)

[~,idx]=sortrows([ { events.year }'  {events.type}' ],[2,1]);
events=events(idx);

% Readjust ANT and REC event to tfirst (take care of duplicates)

for k=1:numel(events)
   if events(k).year < tfirst && any(strcmp(events(k).type,{'ANT','REC'}))
       events(k).year=tfirst;
   end
end

rem=false(numel(events));
for k=2:numel(events)
  if  events(k).year == events(k-1).year && strcmp(events(k).type,events(k-1).type)  
    rem(k-1)=true;
  end
end
events(rem)=[];

% Remove entries before or after selected time

kk=0;
for k=1:numel(events)
  if  events(k).year >= tfirst && events(k).year <= tlast
    kk=kk+1;
    updevents(kk)=events(k);
  end
end

% Set ANT and REC entries to tfirst, if tfirst is earlier than the first
% entry, assuming the early data is the same antenna and receiver type

if tfirst > -Inf
  idxant=find(ismember({ updevents.type },'ANT'));
  if numel(idxant) > 0 && updevents(idxant(1)).year > tfirst
    if opt.verbose > 0
      fprintf('trimEvents: First ANT start time %s adjusted from %s to %s\n',opt.station,datestr(dyear2date(updevents(idxant(1)).year),'yyyy-mm-dd'),datestr(dyear2date(tfirst),'yyyy-mm-dd'));
    end
    updevents(idxant(1)).year=tfirst;
  end
  idxrec=find(ismember({ updevents.type },'REC'));
  if numel(idxrec) > 0 && updevents(idxrec(1)).year > tfirst
    if opt.verbose > 0
       fprintf('trimEvents: First REC start time %s adjusted from %s to %s\n',opt.station,datestr(dyear2date(updevents(idxrec(1)).year),'yyyy-mm-dd'),datestr(dyear2date(tfirst),'yyyy-mm-dd'));
    end
    updevents(idxrec(1)).year=tfirst;
  end
end

end

