function stmdisp(st,varargin)
%STMDISP    Display space time matrix/stuctures/files information.
%  STMDISP(ST) displays information on the space time matrix structure
%  or space time matrix dataset ST. ST can be a space time matrix structure 
%  or the name of a space time matrix dataset. 
%
%  STMDISP(ST,OPTION,VALUE,...) accepts additional OPTION/VALUE pairs.
%  Valid OPTION's are:
%
%    'maxPoints'        % maximum number of points for display
%    'maxEpochs'        % maximum number of epochs for display
%
%  Examples:
%
%    stmdisp('simtest0b_sarAsc1.mat')
%    stmdiff(st)
%
%  See also stm, stmread and stmdiff.
%
%  (c) Hans van der Marel, Delft University of Technology, 2021.

% Created:  14 March 2020 by Hans van der Marel
% Modified: 28 Sep 2023 by Freek van Leijen
%           - allow empty techniqueAttrib in stminfo
%           24 Oct 2023 by Hans van der Marel
%           - added display of atrributes and counts

% Check input arguments and process options

if nargin < 1
   error('This function expects at least one input argument.')
end

opt.maxPoints=50;         % maximum number of points for display
opt.maxEpochs=35;         % maximum number of epochs for display

for k=1:2:numel(varargin)
    if isfield(opt,varargin{k})
       opt.(varargin{k})=varargin{k+1};
    else
       warning(['Invalid option ' varargin{k} ', ignore and continue...'])
    end
end

% If necessary, load space time matrix

if ischar(st)
   stfilename=st;
   st=stmread(stfilename,'NOSTOCHDATA');
   [~,stfilename]=fileparts(stfilename);
elseif isstruct(st)
   stfilename='<struct>';
else
   error('The single input must be a character string with the filename, structure or numeric matrix.')
end

% Print basic info

stminfo(stfilename,st)
printStochModel(st)

% Print point and epoch attributes with counts

stmprinttables(st,opt)

% Print input datasets

printInputDatasets(st.inputDatasets)

% Print global attributes

fprintf('\nGlobal Attributes:\n')
disp(st.globalAttrib)

end

%% Internal functions

function stminfo(stfilename,st)
% Internal function to print basic information

datasetId=st.datasetId;
fprintf('\n%s\n',repmat('-',size(datasetId)))
fprintf('%s\n',datasetId)
fprintf('%s\n\n',repmat('-',size(datasetId)))

datasetAttrib=st.datasetAttrib;

fprintf('datasetId          %s\n', st.datasetId)
fprintf('datasetAttrib:\n')
fprintf('- readAs           %s\n', stfilename)
fprintf('- createdAs        %s\n', datasetAttrib.createdAs)
fprintf('- creationDate     %s\n', datasetAttrib.creationDate)
fprintf('- createdBy        %s\n', datasetAttrib.createdBy)
fprintf('- softwareName     %s\n', datasetAttrib.softwareName)
fprintf('- projectId        %s\n', datasetAttrib.projectId)
fprintf('- projectFile      %s\n', datasetAttrib.projectFile)
fprintf('- projectFileDate  %s\n', datasetAttrib.projectFileDate)
fprintf('\n')

techniqueAttrib=st.techniqueAttrib;

fprintf('techniqueId        %s\n', st.techniqueId)
fprintf('techniqueAttrib:\n')
if ~isempty(techniqueAttrib)
  tmp=fieldnames(techniqueAttrib);
  for k=1:numel(tmp)
     if ischar(techniqueAttrib.(tmp{k})) || isstring(techniqueAttrib.(tmp{k})) 
        fprintf('- %-16s %s\n',tmp{k},techniqueAttrib.(tmp{k}));
     else
        fprintf('- %-16s %s',tmp{k},evalc('disp(st.techniqueAttrib.(tmp{k}))'));
     end
  end
else
  fprintf('-\n');
end
fprintf('\n')

pntAttrib=st.pntAttrib;
maxcrd=max(st.pntCrd);
mincrd=min(st.pntCrd);

fprintf('numPoints          %d\n', st.numPoints)
if ~isempty(pntAttrib)
  fprintf('pntAttrib''s        %s\n', strjoin(fieldnames(pntAttrib),', '));
end
fprintf('latitudeRange      %7.4f - %7.4f\n',mincrd(1),maxcrd(1));
fprintf('longitudeRange     %7.4f - %7.4f\n',mincrd(2),maxcrd(2));
fprintf('\n')

epochAttrib=st.epochAttrib;
yearrange=[ min(st.epochDyear) max(st.epochDyear) ];

fprintf('numEpochs          %d\n', st.numEpochs)
if ~isempty(epochAttrib)
   fprintf('epochAttrib''s      %s\n', strjoin(fieldnames(epochAttrib),', '));
end
fprintf('yearRange          %8.3f - %8.3f\n',yearrange);
fprintf('dateRange          %s - %s\n',datestr(dyear2date(yearrange(1))),datestr(dyear2date(yearrange(2))));
fprintf('\n')

fprintf('numObsTypes        %d\n', numel(st.obsTypes))
fprintf('obsTypes           %s\n', strjoin(st.obsTypes,', '));
fprintf('parTypes           %s\n', strjoin(st.parTypes,', '));
if ~isempty(st.auxTypes)
   fprintf('auxTypes           %s\n', strjoin(st.auxTypes,', '));
end
fprintf('\n')

end

function stmprinttables(st,opt)
% Internal function to print point and epoch tables with counts

if nargin < 2
   opt=struct();
end
if ~isfield(opt,'maxPoints') 
   opt.maxPoints=50;
end
if ~isfield(opt,'maxEpochs') 
   opt.maxEpochs=35;
end
if ~isfield(opt,'pntMask') || isempty(opt.pntMask)
   pntMask=true(size(st.pntName));
end
if ~isfield(opt,'epochMask') || isempty(opt.epochMask)
   epochMask=true(size(st.epochDyear));
end

% Set up point table with counts

pntTable=table();
pntTable.pntName=st.pntName;
pntTable.pntCrd=st.pntCrd;

pntAttribs=fieldnames(st.pntAttrib);
for k=1:numel(pntAttribs)
    pntTable.(pntAttribs{k})=st.pntAttrib.(pntAttribs{k});
end

pntTable.numEpochs=squeeze(sum(~isnan(st.obsData(:,epochMask,:)),2));
pntTable.index=(1:numel(st.pntName))';

pntTable=sortrows(pntTable,'numEpochs','descend');

% Set up point table with counts

if st.numPoints <= opt.maxPoints
   fprintf('\nPoints, attributes and epoch count (sorted on counts):\n\n')
   disp(pntTable)
else
   maxPoints2=ceil(opt.maxPoints/2);
   fprintf('\nPoints, attributes and epoch count (sorted on counts):\n')
   fprintf('- Top %d (out of %d) with highest counts:\n\n',maxPoints2,st.numPoints)
   disp(pntTable(1:maxPoints2,:))
   fprintf('- Bottom %d (out of %d) with lowest counts:\n\n',maxPoints2,st.numPoints)
   disp(pntTable(end-maxPoints2:end,:))
end

% Histogram of epoch counts

if numel(st.obsTypes) > 1
   tmpTable=splitvars(pntTable,'numEpochs');
   histTable=groupcounts(tmpTable,'numEpochs_1');
   histTable=renamevars(histTable,{'numEpochs_1','GroupCount','Percent'}, ...
       {'numEpochs',['pointCount_',st.obsTypes{1}],['Percent_',st.obsTypes{1}]});
   for k=2:numel(st.obsTypes)
      histTable2=groupcounts(tmpTable,['numEpochs_' num2str(k)]);
      histTable2=renamevars(histTable2,{['numEpochs_' num2str(k)'],'GroupCount','Percent'}, ...
       {'numEpochs',['pointCount_',st.obsTypes{k}],['Percent_',st.obsTypes{k}]});
      histTable=outerjoin(histTable,histTable2,'MergeKeys',true);
   end
else
   histTable=groupcounts(pntTable,'numEpochs');
end 
histTable=sortrows(histTable,'numEpochs','descend');

fprintf('\nHistogram counting points with numEpochs observations:\n\n')
disp(histTable)

% Set up epoch table with counts

epochTable=table();
epochTable.epochDyear=st.epochDyear';

epochAttribs=fieldnames(st.epochAttrib);
for k=1:numel(epochAttribs)
    if size(st.epochAttrib.(epochAttribs{k}),1) == st.numEpochs
       epochTable.(epochAttribs{k})=st.epochAttrib.(epochAttribs{k});
    else
      epochTable.(epochAttribs{k})=st.epochAttrib.(epochAttribs{k})';
    end
end

pntMask2 = pntMask & squeeze(sum(~isnan(st.obsData(:,epochMask,:)),2) > 1);

numPnts=zeros(size(st.obsData,2),size(st.obsData,3));
numOk=numPnts;
for l=1:size(numPnts,2)
   numPnts(:,l)=sum(~isnan(st.obsData(pntMask,:,l)),1)';
   numOk(:,l)=sum(~isnan(st.obsData(pntMask2(:,l),:,l)),1)';
end
epochTable.numPnts=numPnts;
epochTable.numOk=numOk;

% Check if epochs are sorted

epochTable.index=(1:numel(st.epochDyear))';
if ~issorted(st.epochDyear)
    disp('ALERT: Epoch data is not sorted on date!!')
    epochTable=sortrows(epochTable,'epochDyear','ascend');
end
epochTable.dDays= [ 0 ; diff(dyear2date(epochTable.epochDyear))];

% Print epoch table

if st.numEpochs <= opt.maxEpochs
   fprintf('\nEpochs, attributes and point count (sorted on time, numOk is numPnts with more than 1 observation):\n\n')
   disp(epochTable)
else
   maxEpochs2=ceil(opt.maxEpochs/2);
   fprintf('\nEpochs, attributes and point count (sorted on time, numOk is numPnts with more than 1 observation):\n')
   fprintf('- Top %d (out of %d) with highest counts:\n\n',maxEpochs2,st.numEpochs)
   disp(epochTable(1:maxEpochs2,:))
   fprintf('- Bottom %d (out of %d) with lowest counts:\n\n',maxEpochs2,st.numEpochs)
   disp(epochTable(end-maxEpochs2:end,:))
end

end

function varargout=printInputDatasets(inputDatasets,indent)
% Internal recursive function to print stm input datasets

if nargin <=1 
    indent=0;
    fprintf('\n')
    fprintf('Input datasets:\n\n')
    fprintf('datasetId                                                                        techniqueId softwareName         creationDate    numPoints numEpochs\n')
    fprintf('-------------------------------------------------------------------------------- ----------- -------------------- --------------- --------- ---------\n')
end

if numel(inputDatasets) > 0
   for k=1:numel(inputDatasets)
      creationDate='';
      if isfield(inputDatasets(k).datasetAttrib,'creationDate')
          creationDate=inputDatasets(k).datasetAttrib.creationDate;
      end
      fprintf('%-80.80s %-11.11s %-20.20s %-15.15s %9d %9d\n', ...
          [ blanks(indent) inputDatasets(k).datasetId  ], ...
          inputDatasets(k).techniqueId, ... 
          inputDatasets(k).datasetAttrib.softwareName, ...
          creationDate, ...
          inputDatasets(k).numPoints, ...
          inputDatasets(k).numEpochs);
      printInputDatasets(inputDatasets(k).inputDatasets,indent+3);
   end
   indent=indent+3;
end

if nargin >1
   varargout={indent};
end

end




