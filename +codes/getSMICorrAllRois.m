% Calculate spatial modulation Z-score Do shuffling on the entire 3D binned
% ROI array
% 
% Input: trial matrix: M(trials,positionBins,ROI)
%           signThreshold default = 0.995;
% Outputs: spatial modulation index (z-score), location of peak 
%
% M = sData.imdata.binnedRoisDeconvRate(sData.trials.contextsMeta(1).trials,:,48);
% M = sData.imdata.binnedRoisDff(sData.trials.contextsMeta(1).trials,:,38);

function [corrCoef, isSignificant, SMI] = getSMICorrAllRois(M1,M2,signThreshold)

if nargin < 3
    signThreshold = 0.995;
end

smoothSpan = 5;

M1 = fillmissing(M1,'linear',2);
MRand1 = M1;
M2 = fillmissing(M2,'linear',2);
MRand2 = M2;


corrCoefShuf = nan(1000,1);

for j = 1:1:1000

% randomize matrix
for r = 1:1:size(M1,1)
    temp = M1(r,:,:);
    i = randi(size(M1,2));
    
    if i ~= 1
        MRand1(r,:,:) = [temp(1,i:size(M1,2),:) temp(1,1:i-1,:)];
    else
        MRand1(r,:,:) = temp;
    end
end
% randomize matrix
for r = 1:1:size(M2,1)
    temp = M2(r,:,:);
    i = randi(size(M2,2));
    
    if i ~= 1
        MRand2(r,:,:) = [temp(1,i:size(M2,2),:) temp(1,1:i-1,:)];
    else
        MRand2(r,:,:) = temp;
    end
end

% even = 3:2:size(M,1); 
% odd = 2:2:size(M,1); 

% Analyze tuning curve

tunCurves1 = smoothdata(permute(nanmean(MRand1),[2 3 1]),'gaussian',smoothSpan);
tunCurves2 = smoothdata(permute(nanmean(MRand2),[2 3 1]),'gaussian',smoothSpan);

for i = 1:1:size(tunCurves1,2)
    corrCoefShuf(j,i) = corr(tunCurves1(:,i),tunCurves2(:,i));
end

end

tunCurve1 = smoothdata(permute(nanmean(M1),[2 3 1]),'gaussian',smoothSpan);
tunCurve2 = smoothdata(permute(nanmean(M2),[2 3 1]),'gaussian',smoothSpan);

corrCoef = []; isSignificant = []; SMI = [];

for i = 1:1:size(tunCurve1,2)
    
    corrCoef(i) = corr(tunCurve1(:,i),tunCurve2(:,i));
    isSignificant(i) = corrCoef(i) > quantile(corrCoefShuf(:,i),signThreshold);
    SMI(i) = (corrCoef(i) - mean(corrCoefShuf(:,i)))/ std(corrCoefShuf(:,i));

end




end


