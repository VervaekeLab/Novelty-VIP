function sData = imQualityChecks(sData)


if ~isfield(sData.imdata, 'roiStat')
    roiStat = codes.getRoiActivityStats(sData);
    sData.imdata.roiStat = roiStat;
end

%% Create roiMeta structure

nROIs = sData.imdata.nROIs;

if ~isfield(sData.imdata,'roiMeta')
    sData.imdata.roiMeta(1:nROIs) = struct();
end

% SNR, baseline...
% sData.imdata(i).roiMeta(roi).SNR = ;

for roi = 1:1:nROIs
    sData.imdata.roiMeta(roi).SNR = sData.imdata.roiStat.signalToNoise(roi);
    sData.imdata.roiMeta(roi).peakDff = sData.imdata.roiStat.peakDff(roi);
    sData.imdata.roiMeta(roi).activityLevel = sData.imdata.roiStat.activityLevel(roi);
end





end