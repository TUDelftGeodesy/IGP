
if isfolder('../stmmain')
   addpath('../stmmain');
   addpath('../stmutil');
elseif isfolder('../0_toolboxes/stmutil')
   addpath('../0_toolboxes/stmutil');
end
options=[];
opt=[];
outputfilename='tt';
flag='overwrite';

% Load space time datasets 
simid='simtest2b';
inputfilenames=[ simid '_filenames_truth.txt'];
[inputfilenames,outputfilename,opt]= ...
    stmcheckarguments(inputfilenames,outputfilename,opt,options,flag);
datasets=[];
for k=1:numel(inputfilenames)
   datasets{k}=stmread(inputfilenames{k});
end

% Select the techniques to process

%techniqueIds=cellfun(@(x) x.techniqueId,datasets,'UniformOutput',false);
%datasetIds=cellfun(@(x) x.datasetId,datasets,'UniformOutput',false);
%selected=ismember(techniqueIds,opt.selectTechnique) & contains(datasetIds,opt.selectDataset);

%datasets(~selected)=[];

numDatasets=numel(datasets);

%% Print overview of datasets, with datasetId, TechniqueId and projectInfo
    
fprintf('\nDatasetId                              TechId createdBy       creationDate    softwareName    projectId       projectFile (projectFileDate)\n')
fprintf(  '-------------------------------------  ------ --------------- --------------- --------------- --------------- -------------------------------\n')
for k=1:numDatasets
   datasetId=datasets{k}.datasetId;
   if length(datasetId) > 38
      datasetId=[datasetId(1:35) '...'];
   end
   datasetAttrib=datasets{k}.datasetAttrib;
   fprintf('%-38s %-5s  %-15s %-15s %-15s %-15s %s (%s)\n',datasetId,datasets{k}.techniqueId, ...
          datasetAttrib.createdBy,datasetAttrib.creationDate,datasetAttrib.softwareName,datasetAttrib.projectId,datasetAttrib.projectFile,datasetAttrib.projectFileDate);
end
fprintf('\n')


%%

for k=1:numel(datasets)
   printStochModel(datasets{k})
end

%%

opt.overrideStochModel=struct();
opt.addStochModel=struct();

%%

% overrideStochModel = struct( 'insar'   , { {'tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)'} } , ...
%                              'sarAsc2' , { {'WN(sigma=1)'} } , ...
%                              'gnss'    , { {'WN(sigma=1)'} {'WN(sigma=1)'} {'WN(sigma=3)'} } )
% addStochModel = struct('insar' , { {'WN(sigma=1)'} } )

opt.overrideStochModel.insar = { {'tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)'} } ;
opt.overrideStochModel.sarAsc2 = { {'WN(sigma=.5)'} } ;
opt.overrideStochModel.gnss =  { {'WN(sigma=1)'} {'WN(sigma=1)'} {'WN(sigma=3)'} } ;
opt.addStochModel.insar = { {'WN(sigma=1)'} };
opt.addStochModel.lev = { {'WN(sigma=1)'} };

%%

for k=1:numel(datasets)
   techniqueId=datasets{k}.techniqueId;
   datasetId=datasets{k}.datasetId;
   stochModel=datasets{k}.stochModel;
   obsTypes=datasets{k}.obsTypes;
   [stochModel, modified] = updStochModel(stochModel,techniqueId,datasetId,obsTypes,opt);
end      

