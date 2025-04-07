function stmdiff(stm1,stm2)
%STMDIFF    Diff (compare) two space time matrix/stuctures/files.
%  STMDIFF(STM1,STM2) compares the space time matrix structures STM1 and 
%  STM2, and prints information on the differences. If possible, the rms 
%  difference between the data is printed. The input arguments STM1 and
%  STM2 can be the name of a space time dataset, or a space time structure,
%  or only the actual space time matrix itself. 
%
%  When STM1 and STM2 are structures or files, this function does it best 
%  to find matching points, epochs and observation types for the comparison. 
%  If necessary, it will use the sensitivity matrix of one of the inputs
%  to match agains the other. This is only possible if the other dataset
%  has the parameter types as observation types. 
%
%  When STM1 and STM2 are not stuctures or files, but numeric matrices, they 
%  must be of the same size: the matrices are compared assuming they refer 
%  to the same points, epochs and observation types. 
%
%  STMDIFF(STM) analyzes the space time matrix or structure STM, it 
%  computes offsets and transforms that minimizes the rms content, and
%  prints the information. 
%
%  Examples:
%
%    stmdiff('simtest0b_sarAsc1.mat','simtest0b_sarAsc1_truth.mat')
%    stmdiff('simtest0b_sarAsc1.mat','simtest0b_truth.mat')
%    stmdiff('simtest0b_combined.mat','simtest0b_truth.mat')
%    stmdiff(stmfile1,stmfile2)
%
%    st1=stmread('simtest0b_combined.mat')
%    st2=stmread('simtest0b_truth.mat')
%    stmdiff(st1,st2)
%
%    stmdiff(st1.obsData,st2.obsData)
%
%  See also stm, stmread and stminterpolate.
%
%  (c) Hans van der Marel, Delft University of Technology, 2020.

% Created:  11 Oct 2020 by Hans van der Marel
% Modified: 28 Oct 2020 by Hans van der Marel
%            - select common observation types
%            - added observation type to print output in compare
%            - read stm from file if necessary
%            - more elaborate statistics on the differences

% Check the input arguments

if nargin < 1
   error('This function requires at least two input arguments.')
end

% Directly call internal function stmstat on the single input

if nargin < 2
   if isstruct(stm1)
      stmstat(stm1.obsData,stm1.obsTypes)
   elseif ischar(stm1)
      stm1=stmread(stm1);
      stmstat(stm1.obsData,stm1.obsTypes)
   elseif isnumeric(stm1)
      stmstat(stm1)
   else
      error('The single input must be a character string with the filename, structure or numeric matrix.')
   end
   return
end

% Directly call internal function stmstat when both input arguments are numeric, and return

if ischar(stm1)
   [~,stmfile1]=fileparts(stm1);
   stm1=stmread(stm1);
else
   stmfile1='STM Dataset 1';
end
if ischar(stm2)
   [~,stmfile2]=fileparts(stm2);
   stm2=stmread(stm2);
else
   stmfile2='STM Dataset 2';
end

if isnumeric(stm1) && isnumeric(stm2)
   if any( size(stm1)-size(stm2) ~= 0 )
       error('Numeric stm matrices must have the same dimensions.')
   end
   stmstat(stm1 - stm2)
   return
elseif ~(isstruct(stm1) && isstruct(stm2))
   error('The inputs must both be structures or both numeric matrices.')
end

% Compare datasetId and techniqueId 

cancompare=true;
comparetype=0;

fprintf('\n                %-30s   %-30s\n',stmfile1,stmfile2)
fprintf('                ------------------------------   ------------------------------\n')

fprintf('datasetId       %-30s | %-30s\n', stm1.datasetId,stm2.datasetId)
fprintf('techniqueId     %-30s | %-30s\n', stm1.techniqueId,stm2.techniqueId)

% Compare the points

if stm1.numPoints == stm2.numPoints
   [d,i1] = pdist2(stm1.pntCrd(:,1:2),stm2.pntCrd(:,1:2),'cityblock','Smallest',1);
   i2=1:size(stm2.pntCrd,1);
   cstr='=';
elseif stm1.numPoints > stm2.numPoints
   [d,i1] = pdist2(stm1.pntCrd(:,1:2),stm2.pntCrd(:,1:2),'cityblock','Smallest',1);
   i2=1:size(stm2.pntCrd,1);
   cstr='>';
elseif stm1.numPoints < stm2.numPoints
   [d,i2] = pdist2(stm2.pntCrd(:,1:2),stm1.pntCrd(:,1:2),'cityblock','Smallest',1);
   i1=1:size(stm1.pntCrd,1);
   cstr='<';
end
maxdiff=max(abs(stm1.pntCrd(i1,:) - stm2.pntCrd(i2,:)));
if all(  maxdiff <1e-5 ) 
   if stm1.numPoints == stm2.numPoints
      fprintf('numPoints       %-30d %1s %-30d  -- Matching coordinates --\n', stm1.numPoints, cstr, stm2.numPoints)
   else
      fprintf('numPoints       %-30d %1s %-30d  -- Subset of %d points with matching coordinates --\n', stm1.numPoints, cstr, stm2.numPoints,numel(d))
   end
else
   fprintf('numPoints       %-30d %1s %-30d  ** Lat, Lon, Hgt diff %.6f %.6f %.2f **\n', stm1.numPoints, cstr, stm2.numPoints, maxdiff)
   cancompare=false;
end

% Compare the epochs

if stm1.numEpochs == stm2.numEpochs
   [d,j1] = pdist2(stm1.epochDyear(:),stm2.epochDyear(:),'cityblock','Smallest',1);
   j2=1:size(stm2.epochDyear(:),1);
   cstr='=';
elseif stm1.numEpochs > stm2.numEpochs
   [d,j1] = pdist2(stm1.epochDyear(:),stm2.epochDyear(:),'cityblock','Smallest',1);
   j2=1:size(stm2.epochDyear(:),1);   
   cstr='>';
elseif stm1.numEpochs < stm2.numEpochs
   [d,j2] = pdist2(stm2.epochDyear(:),stm1.epochDyear(:),'cityblock','Smallest',1);
   j1=1:size(stm1.epochDyear(:),1);   
   cstr='<';
end

maxdiff=max(abs(stm1.epochDyear(j1) - stm2.epochDyear(j2) ));
if all(  maxdiff < 1e-5 ) 
   if stm1.numEpochs == stm2.numEpochs
      fprintf('numEpochs       %-30d %1s %-30d  -- Matching dates --\n', stm1.numEpochs, cstr, stm2.numEpochs)
   else
      fprintf('numEpochs       %-30d %1s %-30d  -- Subset of %d epochs with matching dates --\n', stm1.numEpochs, cstr, stm2.numEpochs,numel(d))
   end
else
   fprintf('numEpochs       %-30d %1s %-30d  ** decimal year diff %.6f **\n', stm1.numEpochs, cstr, stm1.numEpochs, maxdiff)
   cancompare=false;
end

% Compare observation and parameter types

obsTypes1=stm1.obsTypes;
obsTypes2=stm2.obsTypes;

[commonObsTypes,l1,l2]=intersect(obsTypes1,obsTypes2,'stable');

if numel(obsTypes1) == numel(obsTypes2) && numel(obsTypes1) == numel(commonObsTypes)
   cstr='=';
   cmessage='-- Matching observation types --';
elseif numel(commonObsTypes) > 1
   if numel(obsTypes1) < numel(obsTypes2) 
      cstr='<';
   elseif numel(obsTypes1) > numel(obsTypes2) 
      cstr='>';
   else
      cstr='~';
   end
   cmessage=sprintf('-- Subset of %d matching observation types --',numel(commonObsTypes));
else
   cstr='|';
   % check if we can do a comparison in one of the observation domains using the sensitivity matrix
   [parTypesA,l1A,l2A]=intersect(obsTypes1,stm2.parTypes,'stable');
   [parTypesB,l2B,l1B]=intersect(obsTypes2,stm1.parTypes,'stable');
   if numel(parTypesA) >= 3
      comparetype=1;
      cmessage='-- Compare in observation domain of stm2 --';
   elseif numel(parTypesB) >= 3
      comparetype=2;
      cmessage='-- Compare in observation domain of stm1 --';
   else
      cancompare=false;
      cmessage='';
   end
end

fprintf('obsTypes        %-30s %1s %-30s  %s\n', strjoin(stm1.obsTypes),cstr,strjoin(stm2.obsTypes),cmessage)

[commonParTypes]=intersect(stm1.parTypes,stm2.parTypes,'stable');

if numel(stm1.parTypes) == numel(stm2.parTypes) && numel(stm1.parTypes) == numel(commonParTypes)
   cstr='=';
   cmessage='';
elseif numel(commonParTypes) > 1
   cstr='|';
   cmessage=sprintf('-- Subset of %d matching parameter types (this is somewhat unusual) --',numel(commonParTypes));
else
   cstr='|';
   cmessage='** Non matching parameter types (this is unexpected)!! **';
end
fprintf('parTypes        %-30s %1s %-30s  %s\n', strjoin(stm1.parTypes),cstr,strjoin(stm2.parTypes),cmessage)


% Compare the data
%
% - same dimensions (but different order possible) and (some) matching observation types -> compare matching observation types
% - different dimensions but some (matching) observation types -> but subset of points and epochs is close -> select rows and columns
% - non matching observation types (but one is NEU) -> use sensitivity matrix
% - no matching points or epochs  -> use stminterpolate to interpolate (but only when selected as option)
%
% At the moment, only the first three options are implemented

if cancompare
   if comparetype == 0
      stmstat(stm1.obsData(i1,j1,l1) - stm2.obsData(i2,j2,l2),commonObsTypes)
   elseif comparetype == 1
      % [obsTypesA,l1A,l2A]=intersect(obsTypes1,stm2.parTypes,'stable');
      tmp=zeros(numel(i2),numel(j2),numel(obsTypes2));
      for iType=1:numel(obsTypes2)
         for iCol=1:size(stm2.sensitivityMatrix(:,l2A),2)
            d=repmat(stm2.sensitivityMatrix(i2,l2A(iCol),iType),[1 numel(j2)]);
            tmp(:,:,iType)=tmp(:,:,iType)+d.*stm1.obsData(i1,j1,l1A(iCol));   
         end
      end      
      stmstat(tmp - stm2.obsData(i2,j2,:),obsTypes2)
   elseif comparetype == 2
      % [obsTypesB,l2B,l1B]=intersect(obsTypes2,stm1.parTypes,'stable');
      tmp=zeros(numel(i1),numel(j1),numel(obsTypes1));
      for iType=1:numel(obsTypes1)
         for iCol=1:size(stm1.sensitivityMatrix(:,l1B),2)
            d=repmat(stm1.sensitivityMatrix(i1,l1B(iCol),iType),[1 numel(j1)]);
            tmp(:,:,iType)=tmp(:,:,iType)+d.*stm2.obsData(i2,j2,l2B(iCol));   
         end
      end      
      stmstat(stm1.obsData(i1,j1,:) - tmp,obsTypes1)
   else
      error('Hey, this should not happen')
   end
else
   fprintf('\nNo matching points, epochs and/or observation types, cannot compare values.\n')
end

end

function stmstat(stmtmp,obsTypes)

[np,ne,nc]=size(stmtmp);

fprintf('\n')
fprintf('            ___ transf (Epochs) ____      ____ offset (Points) ___      ____________________ difference __________________    (#, eps)\n')
fprintf('            median      min      max      median      min      max      median      mad    stdev      min      max   count\n')
%fprintf('\n           meanEpochs meanPoints   mstdEpochs mstdPoints\n')
for ic=1:nc
   eps=Inf;
   l=0;
   e=zeros(1,ne);
   p=zeros(np,1);
   while eps > 1e-3 && l < 10 
      a=nanmean(stmtmp(:,:,ic),2);
      stmtmp(:,:,ic)=stmtmp(:,:,ic)-repmat(a,[1 size(stmtmp(:,:,ic),2)]);
      b=nanmean(stmtmp(:,:,ic),1);
      stmtmp(:,:,ic)=stmtmp(:,:,ic)-repmat(b,[size(stmtmp(:,:,ic),1) 1]);
      eps=nansum(a.^2)+nansum(b.^2);
      p=p+a;
      e=e+b;
      l=l+1;
   end
   %meanTransf=nanmean(e(:));
   medianTransf=nanmedian(e(:));
   minTransf=nanmin(e(:));
   maxTransf=nanmax(e(:));
   %meanOffset=nanmean(p(:));
   medianOffset=nanmedian(p(:));
   minOffset=nanmin(p(:));
   maxOffset=nanmax(p(:));
   tmp=stmtmp(:,:,ic);
   %meanDispl=nanmean(tmp(:));
   medianDispl=nanmedian(tmp(:));
   stdDispl=nanstd(tmp(:));
   madDispl=1.4826*nanmedian(abs(tmp(:)));
   minDispl=nanmin(tmp(:));
   maxDispl=nanmax(tmp(:));
   cntDispl=sum(~isnan(tmp(:)));
   %mstdEpochs=nanmean(nanstd(stmtmp(:,:,ic),[],1));
   %mstdPoints=nanmean(nanstd(stmtmp(:,:,ic),[],2));
   if nargin < 2
%      fprintf('   %3d    %9.3f  %9.3f    %9.3f  %9.3f    (%d,%9.2e)\n',ic,stdTransf,stdOffsets,mstdEpochs,mstdPoints,l,eps)
      fprintf('   %3d   %9.3f%9.3f%9.3f   %9.3f%9.3f%9.3f   %9.3f%9.3f%9.3f%9.3f%9.3f%8d    (%d,%9.2e)\n',ic, ...
         medianTransf,minTransf,maxTransf, ...
         medianOffset,minOffset,maxOffset, ...
         medianDispl,madDispl,stdDispl,minDispl,maxDispl,cntDispl, ...
         l,eps)
   else
      fprintf('%-9s%9.3f%9.3f%9.3f   %9.3f%9.3f%9.3f   %9.3f%9.3f%9.3f%9.3f%9.3f%8d    (%d,%9.2e)\n',obsTypes{ic}, ...
         medianTransf,minTransf,maxTransf, ...
         medianOffset,minOffset,maxOffset, ...
         medianDispl,madDispl,stdDispl,minDispl,maxDispl,cntDispl, ...
         l,eps)
   end
end

end


   
