function varargout = ClassifyRoisAddSMI(sessionObject, varargin)
%REMAPOPTOPLOTS Summary of this function goes here
%   Detailed explanation goes here


% % % % % % % % % % % % % % % INSTRUCTIONS % % % % % % % % % % % % % % %
% - - - - - - - - - - You can remove this part - - - - - - - - - - - 
% Instructions on how to use this template: 
%   1) If the session method should have parameters, these should be
%      defined in the local function getDefaultParameters at the bottom of
%      this script.
%   2) Scroll down to the custom code block below and write code to do
%   operations on the sessionObjects and it's data.
%   3) Add documentation (summary and explanation) for the session method
%      above. PS: Don't change the function definition (inputs/outputs)
%
%   For examples: Press e on the keyboard while browsing the session
%   methods. (e) should appear after the name in the menu, and when you 
%   select a session method, the m-file will open.


% % % % % % % % % % % % CONFIGURATION CODE BLOCK % % % % % % % % % % % % 
% Create a struct of default parameters (if applicable) and specify one or 
% more attributes (see nansen.session.SessionMethod.setAttributes) for 
% details.
    
    % Get struct of parameters from local function
    params = getDefaultParameters();
    
    % Create a cell array with attribute keywords
    ATTRIBUTES = {'batch', 'queueable'};   

    
% % % % % % % % % % % % % DEFAULT CODE BLOCK % % % % % % % % % % % % % % 
% - - - - - - - - - - Please do not edit this part - - - - - - - - - - - 
    
    % Create a struct with "attributes" using a predefined pattern
    import nansen.session.SessionMethod
    fcnAttributes = SessionMethod.setAttributes(params, ATTRIBUTES{:});
    
    if ~nargin && nargout > 0
        varargout = {fcnAttributes};   return
    end
    
    % Parse name-value pairs from function input and update parameters
    params = utility.parsenvpairs(params, [], varargin);
    
    
% % % % % % % % % % % % % % CUSTOM CODE BLOCK % % % % % % % % % % % % % % 
% Implementation of the method : Add your code here:    
    
for j = 1:1:length(sessionObject)

ind = find(strcmp({sessionObject(j).DataLocation.Name},'sData')); mousePath = fullfile(sessionObject(j).DataLocation(ind).RootPath,sessionObject(j).subjectID);

try
    sDataDir = sessionObject(j).getSessionFolder('sData');
    file = vr.findFileInFolder(sDataDir,sessionObject(j).sessionID,'mat');
    % Select and load sData file
    sessionID = file{1};
    load(fullfile(sDataDir,sessionID));
catch
    % Select and load sData file
    [sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File',mousePath,'MultiSelect','off' );
    load(fullfile(sDataDir,sessionID));
end

sessionID = strsplit(sessionID,'.mat');
sessionID = sessionID{1};

%% Bootstrap


%     SMI = []; SMIcorr = [];
%     for t = 1:1:length(sData.trials.contextsMeta)
%         trials = sData.trials.contextsMeta(t).trials;
%         dataM = sData.imdata.binnedRoisDff(trials,:,:);
%         
%         [~, ~, SMIcorr(:,t)] = vr.getSMICorrAllRois(dataM(1:2:numel(trials),:,:),dataM(2:2:numel(trials),:,:),0.95);
%        
%         nROIs = sData.imdata.nROIs;
% 
%         for r = 1:1:nROIs
%             SMI(r,t) = vr.getSMI(dataM(:,:,r));
%         end
%       
%     end
%     
%     for r = 1:1:nROIs
%         
%         sData.imdata.roiMeta(r).smiPeak = SMI(r,:);
%         sData.imdata.roiMeta(r).smiCorr = SMIcorr(r,:);
%         
%     end
 
%% Find place fields

% disp('Place cells:' )
% disp([sData.imdata.placeCells.nPlaceCells])
% disp('Place cell fractions:' )
% disp([sData.imdata.placeCells.placeCellFraction])

% % % % % disp('Before:')
% % % % % disp([[sData.imdata.placeCells.nPlaceCells], round([sData.imdata.placeCells.placeCellFraction]*100)])
% % % % % 

%sData = vr.iman.classifyROIsV3(sData,sDataDir);

% % % % % 
% % % % % disp('After:')
% % % % % disp([[sData.imdata.placeCells.nPlaceCells], round([sData.imdata.placeCells.placeCellFraction]*100)])



%%
%% Control Figure

Xax = sData.stats.sessionAvs(1).plotXAxis;

fov = 1;

figure('Color','white','Position',[0 0 1200 800])

placeCells = reshape([sData.imdata.roiMeta(:).placeCell], length(sData.trials.contextsMeta),[]); 
for b = 1:1:length(sData.trials.contextsMeta)
    
trialsOdd = sData.trials.contextsMeta(b).trials(1:2:sData.trials.contextsMeta(b).nTrials);
trialsEven = sData.trials.contextsMeta(b).trials(2:2:sData.trials.contextsMeta(b).nTrials);

posTunCurvesOdd =  nanmean(sData.imdata(fov).binnedRoisDff(trialsOdd,:,:)); 
posTunCurvesOdd = permute(posTunCurvesOdd,[3 2 1]);

posTunCurvesEven =  nanmean(sData.imdata(fov).binnedRoisDff(trialsEven,:,:)); 
posTunCurvesEven = permute(posTunCurvesEven,[3 2 1]);
posTunCurvesEven =  smoothdata(posTunCurvesEven,2,'gaussian',5);

%posTunCurves = sData.imdata(fov).avBinnedRois.avBinnedRoisDff{b};
%posTunCurves =  smoothdata(posTunCurves,2,'gaussian',smoothSpan);
nROIs = sData.imdata.nROIs;

cMin = quantile(posTunCurvesEven(:),0.01);
cMax = quantile(posTunCurvesEven(:),0.99);
cLabel = 'Dff';





sorted = vr.sortROIs(posTunCurvesOdd,find(placeCells(b,:)),5);

subplot(2,length(sData.trials.contextsMeta),b)
%subplot(1,2,1)
imagesc(Xax,1:numel(sorted),posTunCurvesEven(sorted,:))
%title('Place cells') 
t = title('Place cells');
t.FontWeight = 'normal';
t.FontSize = 14;
%c = colorbar; c.Label.String = cLabel;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
%xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
caxis([cMin cMax]);


sorted = vr.sortROIs(posTunCurvesOdd,setdiff(1:nROIs,find(placeCells(b,:))),5);
subplot(2,length(sData.trials.contextsMeta),length(sData.trials.contextsMeta) + b)
%subplot(1,2,2)
imagesc(Xax,1:numel(sorted),posTunCurvesEven(sorted,:))
%title('Non-place cells') 
t = title('Non-place cells');
t.FontWeight = 'normal';
t.FontSize = 14;
%c = colorbar; c.Label.String = cLabel;
colormap(gca,jet);
ax = gca;
ax.TickDir = 'out';
xlabel('\fontsize{13}Track position (cm)');
y = ylabel('\fontsize{13}sorted ROIs');
%y.FontWeight = 'bold';
caxis([cMin cMax]);


end
suptitle([sData.sessionInfo.sessionID(1:17) ' - ' 'Fov' num2str(fov) ' ' sData.imdata(fov).fovLocation ' - ' sData.trials.vrContextProtocols(b).name])

disp('Place cells:' )
disp([sData.imdata.placeCells.nPlaceCells])
disp('Place cell fractions:' )
disp([sData.imdata.placeCells.placeCellFraction])






%% Save sData

save(fullfile(sDataDir,sessionID),'sData');

disp([sessionID ' - ROI classification done'])





end
    
    
    % Return session object (please do not remove):
    % if nargout; varargout = {sessionObject}; end
end


function params = getDefaultParameters()
%getDefaultParameters Get the default parameters for this session method
%
%   params = getDefaultParameters() should return a struct, params, which 
%   contains fields and values for parameters of this session method.

    % Add fields to this struct in order to define parameters for this
    % session method:
    params = struct();

end