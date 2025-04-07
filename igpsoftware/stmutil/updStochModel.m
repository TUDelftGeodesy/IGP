function [stochModel, modified] = updStochModel(stochModel,techniqueId,datasetId,obsTypes,opt)
%updStochModel  Update stochastic model for space time datasets.
%   stochModel = updStochModel(stochModel,techniqueId,datesetId,numObsTypes, opt)
%   updates the stochastic model parameters for space time datasets. The 
%   input arguments are the stochModel, techniqueId, datasetId and numObsTypes 
%   from the space time dataset, and a stucture with replacement data.
%   The following fields of opt are used (if present)
%
%      opt.overrideStochModel   structure with replacement stochastic models
%      opt.addStochModel        structure with stochastic models to add
%      opt.verbose              scalar with verbosity level
%
%   Examples of the overrideStochModel and addStochModel structures are
%
%      opt.overrideStochModel.insar = { {'tudinsar4(s20=7.93,s2t=5.5,s2s=3.9,Rt=0.67,Rs=1.11)'} } ;
%      opt.overrideStochModel.sarAsc2 = { {'WN(sigma=1)'} } ;
%      opt.overrideStochModel.gnss =  { {'WN(sigma=1)'} {'WN(sigma=1)'} {'WN(sigma=3)'} } ;
%
%      opt.addStochModel.insar = { {'WN(sigma=1)'} };
%
%   where the fieldname of the structure must be equal to a valid techniqueId
%   or datasetId. The input arguments techniqueId and datasetId are matched
%   against the structure fieldnames. The values are a cell array with
%   numObsTypes element, with each element a cell (array) with the
%   stochastic model parameters. The values of overrideStochModel and
%   addStochModel are the same format as the values of the input and
%   output stochModel.
%
%   As indicated by the name overrideStochModel replaces the original
%   model(s) in the input stochModel, while addStochModel add's an
%   additional model to the existing stochModel.
%   
%   The output consists of the updated stochModel.
%
%   [stochModel, modified] = stochModelUpd(...) returns also a boolean
%   with true if the model was updated, and false if not.
%
%   See also stm, stmstochmodel and prtStochModel
%
%   (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:  22 December 2020 by Hans van der Marel
% Modified: 27 October 2022 by Hans van der Marel
%              - replace decimal points in datasetId by underscores
%              - replaced ismember by strncmpi in datasetId test and
%                compare only on the shortest string

% Check input arguments 

if nargin ~= 5
   error('This function expects exactly 5 input arguments.')
end

% Replace decimal point in datasetId by underscore

datasetId = strrep(datasetId,'.','_');

% Prepare options and outputs

if isfield(opt,'verbose')
   verbose=opt.verbose;
else
   verbose=1;
end
modified=false;
numObsTypes = numel(obsTypes);

% Optionally override stochModel parameters 

if isfield(opt,'overrideStochModel')

    overrideStochModel=opt.overrideStochModel;
    
    % Override stochModel parameters based on techniqueId
    
    override=ismember(techniqueId,fieldnames(overrideStochModel));
    if any(override)
      if (verbose > 0), fprintf('Override %s stochastic model for dataset %s\n',techniqueId,datasetId); end
      stochModel=overrideStochModel.(techniqueId);
      modified=true;
    end

    % Override stochModel parameters based on datasetId

%     override=ismember(datasetId,fieldnames(overrideStochModel));
%     if any(override) && ~strcmp(datasetId,techniqueId)
%       if (verbose > 0), fprintf('Override %s stochastic model\n',datasetId); end
%       stochModel=overrideStochModel.(datasetId);
%       modified=true;
%     end

    fields=fieldnames(overrideStochModel);
    for k=1:numel(fields)
        if strncmpi(datasetId,fields{k},min(length(datasetId),length(fields{k}))) && ~strcmp(datasetId,techniqueId)
           if (verbose > 0), fprintf('Override %s stochastic model\n',datasetId); end
           stochModel=overrideStochModel.(fields{k});
           modified=true;
        end
    end

end

% Optionally add stochModel parameters 

if isfield(opt,'addStochModel')

    addStochModel=opt.addStochModel;

    % Add model to stochModel based on techniqueId

    add=ismember(techniqueId,fieldnames(addStochModel));
    if any(add)
      if (verbose > 0), fprintf('Add %s stochastic model for dataset %s\n',techniqueId,datasetId); end
      tmp = addStochModel.(techniqueId);
      for l=1:numObsTypes
         stochModel{l}= [ cellstr(stochModel{l}) cellstr(tmp{l}) ];
      end
      modified=true;
    end

    % Add model to stochModel based on datasetId

%     add=ismember(datasetId,fieldnames(addStochModel));
%     if any(add) && ~strcmp(datasetId,techniqueId)
%       if (verbose > 0), fprintf('Add %s stochastic model\n',datasetId); end
%       tmp = addStochModel.(datasetId);
%       for l=1:numObsTypes
%          stochModel{l}= [ cellstr(stochModel{l})  cellstr(tmp{l}) ];
%       end      
%       modified=true;
%     end

    fields=fieldnames(addStochModel);
    for k=1:numel(fields)
        if strncmpi(datasetId,fields{k},min(length(datasetId),length(fields{k}))) && ~strcmp(datasetId,techniqueId)
           if (verbose > 0), fprintf('Override %s stochastic model\n',datasetId); end
           if (verbose > 0), fprintf('Add %s stochastic model\n',datasetId); end
           tmp = addStochModel.(fields{k});
           for l=1:numObsTypes
              stochModel{l}= [ cellstr(stochModel{l})  cellstr(tmp{l}) ];
           end      
           modified=true;
        end
    end
 
end

if modified && verbose > 0
   printStochModel(stochModel,techniqueId,datasetId,obsTypes)
end

end
