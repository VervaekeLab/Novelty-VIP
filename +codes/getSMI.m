% Calculate spatial modulation Z-score 
% Input: trial matrix: M(trials,positionBins)
% Outputs: spatial modulation index (z-score), location of peak 
%
% M = sData.imdata.binnedRoisDeconvRate(sData.trials.trialTypesMeta(1).trials,31:165,38);
% M = binnedRoisDff(sData.trials.trialTypesMeta(2).hitTrials,31:165,15);

function [SMI, peakAmplRand] = getSMI(M)
%???? What should I do with the NaNs?
M = fillmissing(M,'linear',2);
MRand = M;

smoothSpan = 9;
peakAmplRand = NaN(1000,1);
tunCurves = NaN(1000,size(M,2));

for j = 1:1:1000

% randomize matrix
for r = 1:1:size(M,1)
    temp = M(r,:);
    i = randi(size(M,2));
    
    if i ~= 1
        MRand(r,:) = [temp(i:size(M,2)) temp(1:i-1)];
    else
        MRand(r,:) = temp;
    end
end


% Analyze tuning curve

tunCurve = smoothdata(nanmean(MRand),2,'gaussian',smoothSpan);
tunCurves(j,:) = tunCurve;
peakAmplRand(j) = (max(tunCurve) - min(tunCurve)); % / mean(tunCurve);

end

tunCurve = smoothdata(nanmean(M),2,'gaussian',smoothSpan);
peakAmplExp = (max(tunCurve) - min(tunCurve)); % / mean(tunCurve);

peakBin = find(tunCurve == max(tunCurve));

SMI = abs(peakAmplExp-nanmean(peakAmplRand))/std(peakAmplRand);


end

