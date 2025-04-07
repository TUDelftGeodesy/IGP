function printStochModel(stochModel,techniqueId,datasetId,obsTypes)
%printStochModel  Print stochastic model for space time datasets.
%   printStochModel(stochModel,techniqueId,datesetId,numObsTypes) prints
%   the stochastic model parameters for space time datasets. The 
%   input arguments are the stochModel, techniqueId, datasetId and obsTypes 
%   from the space time dataset.
%
%   printStochModel(stm) prints the stochastic model parameters for space
%   time dataset stm, with stm the space time matrix structure.
%
%   See also stm and stmstochmodel.
%
%   (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:  22 December 2020 by Hans van der Marel
% Modified: 

% Check input arguments 

if nargin ~= 1 && nargin ~= 4 
   error('This function expects 1 or 4 input arguments.')
end
if nargin ==1
   % Input is a dataset
   techniqueId=stochModel.techniqueId;
   datasetId=stochModel.datasetId;
   obsTypes=stochModel.obsTypes;
   stochModel=stochModel.stochModel;
end

% Print

numObsTypes=numel(obsTypes);
if numObsTypes > 1
   fprintf('Stochastic model for %s dataset %s\n',techniqueId,datasetId);
else
   fprintf('Stochastic model for %s (%s) dataset %s\n',techniqueId,obsTypes{1},datasetId);
end

for l=1:numObsTypes
   if numObsTypes > 1, fprintf('%s:\n',obsTypes{l}); end
   if l <= numel(stochModel)
      tmpModel=cellstr(stochModel{l});
   else
      tmpModel={'No stochastic model available for this obsType.'};
   end
   for k=1:numel(tmpModel)
      fprintf('   %s\n',tmpModel{k});
   end
end

end