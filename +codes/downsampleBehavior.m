 function sData = downsampleBehavior(sData)

try
    
corridorPosition = sData.behavior.signals.corridorPosition; 
trialStartIndexes = find(diff(corridorPosition) < -10)+1;
nAllTrials = numel(trialStartIndexes)-1;

lickPositions = sData.behavior.signals.lickPositions;

fs = sData.daqdata.meta.fs;

positionRelativeToStart = nan(size(lickPositions));
positionRelativeToFirstLick = nan(size(lickPositions));
timeRelativeToStart = nan(size(lickPositions));
timeRelativeToLandmark = nan(size(lickPositions));
timeRelativeToFirstLick = nan(size(lickPositions));


for i = 1:1:nAllTrials
    
    j = trialStartIndexes(i);
    k = trialStartIndexes(i+1)-1;
    
    tempPos = corridorPosition(j:k);
    tempLickPos = lickPositions(j:k);
    temp = tempLickPos(tempLickPos>0);
        
    tempTime = linspace(0,numel(tempPos)/fs,numel(tempPos));
    landmarkIndex = find(abs(tempPos) == min(abs(tempPos)));
    landmarkIndex = landmarkIndex(1);
    timeRelativeToStart(j:k) = tempTime;
    timeRelativeToLandmark(j:k) = tempTime - landmarkIndex/fs;
    positionRelativeToStart(j:k) = tempPos - tempPos(1);
    
    if numel(temp) > 0
        positionRelativeToFirstLick(j:k) = tempPos - temp(1);
        tempFirstLickIndex = find(tempLickPos == temp(1));
        timeRelativeToFirstLick(j:k) = tempTime - (tempFirstLickIndex(1)/fs);
    end
    
    clear('tempPos', 'tempLickPos', 'temp', 'tempTime', 'landmarkIndex');
end
    
catch
end


%% 

corridorPosition = sData.behavior.signals.corridorPosition; 

trialStartIndexes = find(diff(corridorPosition) < -10) +1; % changed from -100 to capture very first trials
trialIndexes = nan(size(sData.behavior.signals.corridorPosition));

fs = sData.daqdata.meta.fs;


for t = 1:1:numel(trialStartIndexes)-1
    trialIndexes(trialStartIndexes(t):trialStartIndexes(t+1)-1) = t;   
end

sData.behavior.signals.trialIndexes = trialIndexes;

%corridorPosition = sData.behavior.signals.corridorPosition; 




    scanFrameRate = sData.imdata.meta.fps;
    samplePerFrame = fs/scanFrameRate;
    frameIndexes = codes.fixFrameIndexes(sData.daqdata.frameIndex, samplePerFrame);


if numel(frameIndexes) > sData.imdata.nSamples
    frameIndexes = frameIndexes(1:sData.imdata.nSamples);
end


timeDs = frameIndexes/fs;


sData.behavior.signals.trialIndexesDs = trialIndexes(frameIndexes);
sData.behavior.signals.corridorPositionDs = corridorPosition(frameIndexes);

sData.behavior.signals.velocityDs = sData.behavior.signals.velocity(frameIndexes);
sData.behavior.signals.lickIfreqDs = sData.behavior.signals.lickIfreq(frameIndexes);
try
sData.behavior.signals.stimProtsDs = sData.daqdata.stimProtIndicators(frameIndexes); 
catch
end
sData.behavior.signals.timeDs = timeDs; %0:deltaTime:deltaTime*(numel(frameIndexes)-1); %the most accutate way is to calculate from the fram signals


% try
% 
% sData.behavior.signals.timeRelativeToStartDs = timeRelativeToStart(frameIndexes);
% sData.behavior.signals.timeRelativeToLandmarkDs = timeRelativeToLandmark(frameIndexes);
% sData.behavior.signals.timeRelativeToFirstLickDs = timeRelativeToFirstLick(frameIndexes);
% sData.behavior.signals.positionRelativeToStartDs = positionRelativeToStart(frameIndexes);
% sData.behavior.signals.positionRelativeToFirstLickDs = positionRelativeToFirstLick(frameIndexes);
% sData.behavior.signals.rewardDs = sData.daqdata.waterValve(frameIndexes);
% catch
% end


if isfield(sData.daqdata,'optoSignal')
   sData.behavior.signals.optoSignalDs = sData.daqdata.optoSignal(frameIndexes);    
end
    




    
end