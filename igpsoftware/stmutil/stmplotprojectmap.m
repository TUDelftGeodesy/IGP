function stmplotprojectmap(stmProject,stmIntegrate)
%stmplotprojectmap   Plot network of evaluation points.
%   STMPLOTPROJECTMAP(STMPROJECT) plot the network of evaluation points from
%   the project file STMPROJECT showing possible evaluation points and 
%   their origin.
%
%   STMPLOTPROJECTMAP(STMPROJECT, STMINTEGRATE) does the same, but only shows
%   the evaluation points used durring the integration. Needs the integration
%   STM matrix STMINTEGRATE as second input.
%
%   Example:
%
%      stmplotprojectmap('../2_reduce/waddenzee.mat', ...
%             '../3_integrate/waddenzee_insarReducedWithDiag_LevDiag.mat');
%
%   See also STMSELECT and STMINTEGRATE.
%
%  (c) Hans van der Marel, Delft University of Technology, 2023.

% Created:  12 April 2023 by Hans van der Marel
% Modified: 23 Oct 2023 by Hans van der Marel
%            - renamed stmnetworkplot to stmplotprojectmap
%            - check if first stm is a project file
%            - extended datatips

%% Initialize

if nargin < 1 || nargin > 2
    error('stmplotprojectmap expects one or two arguments')
end 
if nargin < 2
    stmIntegrate = '';
end

%% Load project space time matrix

stp=stmread(stmProject);

if ~strcmpi(stp.techniqueId,'projectFile')
    error('First space-time matrix must be a projectFile.')
end

pntCrd = stp.pntCrd;
pntName = stp.pntName;
pntId = stp.pntAttrib.pntId;

markerDatasetIds = stp.techniqueAttrib.markerDatasetIds;
campaignDatasetIds = stp.techniqueAttrib.campaignDatasetIds;

hasMarker = stp.pntAttrib.hasMarker;
dsPntIndex = stp.pntAttrib.dsPntIndex;
dsMarkerNames = stp.pntAttrib.dsMarkerNames;

densificationROI = stp.datasetAttrib.softwareOptions.densificationROI;
inDensificationROI = inpolygon(pntCrd(:,1),pntCrd(:,2),densificationROI(:,1),densificationROI(:,2));


%% Print summary 

fprintf('\nNumber of candidate evaluation points:   %d\n',size(pntCrd,1))
fprintf('- original   %5d\n',sum((hasMarker==1)))
fprintf('- merged     %5d\n',sum((hasMarker>1)))
fprintf('- new        %5d\n',sum((hasMarker==0)))
%epochIntervals=diff(epochDyearTotal);
%fprintf('Min/mean/median/max interval between epochs: %.2f / %.2f / %.2f / %.2f [yr]\n',min(epochIntervals),mean(epochIntervals),median(epochIntervals),max(epochIntervals))
%fprintf('MAD/StDev of intervals:  %.2f / %.2f [yr]\n',1.4826*median(abs( epochIntervals - median(epochIntervals))),std(epochIntervals))

%% Read info from integration stm

if ~isempty(stmIntegrate)

   sti = stmread(stmIntegrate,'NODATA');

   integrationROI = sti.datasetAttrib.softwareOptions.ROI;
   inIntegrationROI = inpolygon(pntCrd(:,1),pntCrd(:,2),integrationROI(:,1),integrationROI(:,2));

   [isUsed, loci] = ismember(stp.pntAttrib.pntId, sti.pntAttrib.pntId);

   fprintf('\nNumber of points (used):   %d\n',sum(isUsed))
   fprintf('- original   %5d\n',sum((hasMarker==1 & isUsed)))
   fprintf('- merged     %5d\n',sum((hasMarker>1 & isUsed)))
   fprintf('- new        %5d\n',sum((hasMarker==0 & isUsed)))

   inROI = inIntegrationROI | inDensificationROI;

   pntName(isUsed)=sti.pntName(loci(isUsed));

else

   integrationROI = densificationROI;
   inROI = inDensificationROI;
   isUsed = inROI;

end

%% Prepare network plot

% Check if datatips are supported

opt.datatips=true;
if opt.datatips && exist('dataTipTextRow','file') ~= 2
    fprintf('Warning: Datatips are not supported (by this Matlab version)')
    opt.datatips=false;
end

% Map labels and aspect ratio

if isfield(opt,'doplots_mapcrd') && opt.doplots_mapcrd
   xylabel = { 'East [km]', 'North [km]' };
   xyaspect = [1 1 1];
else
   xylabel = { 'Longitude [deg]', 'Latitude [deg]' };
   xyaspect = [ 1/cosd(mean(pntCrd(:,1))) 1 1];
end

%% Create network plot

figure;hold on;
for k=1:numel(markerDatasetIds)
    indataset = ~isnan(dsPntIndex(:,k));
    technique = stp.inputDatasets(ismember(  { stp.inputDatasets.datasetId }, markerDatasetIds{k})).techniqueId;
    if ~ismember(markerDatasetIds{k},campaignDatasetIds)
       hs1=plot(pntCrd(indataset & inROI,2),pntCrd(indataset & inROI,1),'^','Markersize',4,'DisplayName',markerDatasetIds{k});
    elseif strcmpi(technique,'gnss')
       hs1=plot(pntCrd(indataset & inROI,2),pntCrd(indataset & inROI,1),'d','Markersize',3,'DisplayName',markerDatasetIds{k});
    else
       hs1=plot(pntCrd(indataset & inROI ,2),pntCrd(indataset & inROI ,1),'.','Markersize',6,'DisplayName',markerDatasetIds{k});
    end
    if opt.datatips
        hs1.DataTipTemplate.DataTipRows = [ ...
            dataTipTextRow('pntName',dsMarkerNames(indataset & inROI,k)), ...
            dataTipTextRow('pntId',pntId(indataset & inROI)), ...
            dataTipTextRow('Lon [deg]',pntCrd(indataset & inROI,2)) , ...
            dataTipTextRow('Lat [deg]',pntCrd(indataset & inROI,1)) ];
    end
end   
hs2=plot(pntCrd(hasMarker>1 & isUsed,2),pntCrd(hasMarker>1 & isUsed,1),'ro','markersize',8,'DisplayName','Merged points');
if opt.datatips
    hs2.DataTipTemplate.DataTipRows = [ ...
       dataTipTextRow('pntName',pntName(hasMarker>1 & isUsed)), ...
       dataTipTextRow('Lon [deg]',pntCrd(hasMarker>1 & isUsed,2)) , ...
       dataTipTextRow('Lat [deg]',pntCrd(hasMarker>1 & isUsed,1)) ];
end
hs3=plot(pntCrd(hasMarker==0 & isUsed,2),pntCrd(hasMarker==0 & isUsed,1),'g+','markersize',3,'DisplayName','InSAR only');
if opt.datatips
    hs3.DataTipTemplate.DataTipRows = [ ...
       dataTipTextRow('pntName',pntName(hasMarker==0 & isUsed)), ...
       dataTipTextRow('Lon [deg]',pntCrd(hasMarker==0 & isUsed,2)) , ...
       dataTipTextRow('Lat [deg]',pntCrd(hasMarker==0 & isUsed,1)) ];
end
plot(densificationROI(:,2),densificationROI(:,1),'-','HandleVisibility', 'off')
plot(integrationROI(:,2),integrationROI(:,1),'-','HandleVisibility', 'off')
xlabel(xylabel{1});
ylabel(xylabel{2});
legend('Interpreter','None','Location','best');

title('Evaluation points (closeby points merged)');

daspect(xyaspect);

end
