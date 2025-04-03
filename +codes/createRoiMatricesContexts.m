% Create trial matrices for all ROIs. Part of image analysis pipeline.

% 
function sData = createRoiMatricesContexts(sData)

% SET VELOCITY THRESHOLD
velThreshold = 0.5; % Do not consider calcium signal if the mouse stops, or walks slower then this value.

%% Prepare data and create trial matrices for all ROIs

samplingRate = sData.daqdata.meta.samplingRate;
binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;
try
    homeBoxLength = abs(sData.stats.sessionAvs(1).plotXAxis(1) - binSize);
catch
    homeBoxLength = 0;
end
%rewardZone = sData.behavior.rewardZone;
corridorLength = binNumber*binSize;
%viewDistance = sData.behavior.viewDistance;

%nHomeBoxBins = homeBoxLength/binSize;
binnedPosition = discretize(sData.behavior.signals(1).corridorPosition,-homeBoxLength:binSize:-homeBoxLength+corridorLength);
trialStartIndexes = find(diff(sData.behavior.signals(1).corridorPosition)< -10)+1;
allTrials = numel(trialStartIndexes)-1;

if allTrials ~= sData.behavior.trialMatrices.meta.nAllTrials % in rare cases the fillmissing smoothes the transition between the trials thererore the raw position data can be better
    trialStartIndexes = find(diff(sData.daqdata.unityPosition)< -10)+1;
    allTrials = numel(trialStartIndexes)-1;
end

%velocity = sData.behavior.signals.velocity;
velocity = smoothdata(sData.behavior.signals(1).velocity,'gaussian',samplingRate/30); % Gaussian filter 30 Hz


%%

% ROI signals 
dff = sData.imdata.roiSignals(2).dff;

% try
% ciaDeconv = sData.imdata.roiSignals(2).deconv;
% deconvDataExist = true;
% catch
%     ciaDeconv = dff;
%     deconvDataExist = false;
% end
% try
% deconvRate = sData.imdata.roiSignals(2).actRateDeconv;
% deconvRateDataExist = true;
% catch
%     deconvRate = ciaDeconv;
%     deconvRateDataExist = false;
% end

nROIs = numel(dff(:,1));


    scanFrameRate = sData.imdata.meta.fps;
    samplePerFrame = sData.daqdata.meta.fs/scanFrameRate;
    frameIndexes = codes.fixFrameIndexes(sData.daqdata.frameIndex, samplePerFrame);


frameSignalInd = zeros(1,numel(sData.daqdata.frameIndex)); % This modified frame signal contains the index in the frame array ..001..0002.. 003

for ind = 1:1:min(numel(frameIndexes),sData.imdata.nSamples)

    frameSignalInd(frameIndexes(ind)) = ind;
    
end


binnedRoisDff = nan(allTrials,binNumber,nROIs);
% binnedRoisDeconv = nan(allTrials,binNumber,nROIs);
% binnedRoisDeconvRate = nan(allTrials,binNumber,nROIs);


% Create roi matrices DFF and Deconv
for r = 1:1:allTrials
    tempPos = binnedPosition(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    tempFrameSInd = frameSignalInd(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    tempVel = velocity(trialStartIndexes(r):(trialStartIndexes(r+1)-1));
    tempFrameSInd(find(tempVel < velThreshold)) = 0;
    tempBinPos = tempPos(find(tempFrameSInd > 1));
    
    for roi = 1:1:nROIs
        tempSignalDff = dff(roi,tempFrameSInd(tempFrameSInd>0));
%        tempSignalDeconv = ciaDeconv(roi,tempFrameSInd(tempFrameSInd>0));
%        tempSignalDeconvRate = deconvRate(roi,tempFrameSInd(tempFrameSInd>0));
        for c = 1:1:binNumber
            binnedRoisDff(r,c,roi) = mean(tempSignalDff(find(tempBinPos == c)));
%            binnedRoisDeconv(r,c,roi) = mean(tempSignalDeconv(find(tempBinPos == c)));  
%            binnedRoisDeconvRate(r,c,roi) = mean(tempSignalDeconvRate(find(tempBinPos == c)));  
        end
        clear('tempSignalDff','tempSignalDeconv','tempSignalDeconvRate');
    end
    clear('tempPos','tempFrameSInd','tempBinPos','tempVel');
end



sData.imdata.binnedRoisDff = fillmissing(binnedRoisDff,'linear',2);
% sData.imdata.binnedRoisDeconv = fillmissing(binnedRoisDeconv,'constant',0);
% if deconvRateDataExist; sData.imdata.binnedRoisDeconvRate = fillmissing(binnedRoisDeconvRate,'linear',2); end



%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Create position tuning curves

binNumber = sData.behavior.trialMatrices.meta.binNumber;
binSize = sData.behavior.trialMatrices.meta.binSize;


nCtxBlocks = length(sData.trials.contextsMeta);

avBinnedRoisDff = {};
%avBinnedRoisDeconv = {};
%avBinnedRoisDeconvRate = {};

for c = 1:1:nCtxBlocks
%avBinnedRoisTrialType1Dff(nROIs,binNumber) = NaN;
for roi = 1:1:nROIs 
    avBinnedRoisDff{c}(roi,:) = nanmean(binnedRoisDff(sData.trials.contextsMeta(c).trials,:,roi));
%    avBinnedRoisDeconv{c}(roi,:) = nanmean(binnedRoisDeconv(sData.trials.contextsMeta(c).trials,:,roi));
%    avBinnedRoisDeconvRate{c}(roi,:) = nanmean(binnedRoisDeconvRate(sData.trials.contextsMeta(c).trials,:,roi));
end
end



sData.imdata.binNumber = binNumber;
sData.imdata.binSize = binSize;

sData.imdata.avBinnedRois.avBinnedRoisDff = avBinnedRoisDff;
%if deconvDataExist; sData.imdata.avBinnedRois.avBinnedRoisDeconv = avBinnedRoisDeconv; end
%if deconvRateDataExist; sData.imdata.avBinnedRois.avBinnedRoisDeconvRate = avBinnedRoisDeconvRate; end


end


% 
% function frameIndexes = fixFrameIndexes(frameSignal, samplePerFrame)
% % This function does a quality check on the frame signal. 
% % In Scanimage (mesoscope recordings) the frame signal occasionally does not 
% % return to 0 which can result in skipping frames. 
% % (check the time stamp in the tiff file!!) 
% % Example session with error: 'm7128-20230322-01-NL4-multiRSC' 
% 
% 
% %Digital
% % frameSignal = sData.daqdata.frameIndex;
% % frameSignal = sData.daqdata.frameIndex2;
% % samplePerFrame = 165;
% 
% if nargin < 2
%     samplePerFrame = median(diff(find(diff(frameSignal)==1)));
% end
% 
% frameIndexes = find(diff(frameSignal)==1);
% 
% incrementSample = diff(frameIndexes);
% 
% incrementFrame = round(incrementSample/samplePerFrame);
% 
% errors =  find(incrementFrame>1);
% 
% for i = flip(1:1:numel(errors))
% 
%     j = errors(i);
%     k = incrementFrame(j);
%     insertion = 1:1:k-1; 
%     insertion = insertion * round(incrementSample(j)/incrementFrame(j)) +frameIndexes(j);
%     frameIndexes = [frameIndexes(1:j), insertion, frameIndexes(j+1:end)];
% 
% end
% 
% if numel(errors) > 0
%     msgbox([num2str(numel(errors)) ' missing frames have been corrected.'])
%     figure
%     histogram(incrementFrame)
%     title('nFrame/frameIndex')
% 
% end
% 
% end