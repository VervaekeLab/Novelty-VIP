function frameIndexes = fixFrameIndexes(frameSignal, samplePerFrame)
% This function does a quality check on the frame signal. 
% In Scanimage (mesoscope recordings) the frame signal occasionally does not 
% return to 0 which can result in skipping frames. 


if nargin < 2
    samplePerFrame = median(diff(find(diff(frameSignal)==1)));
end

frameIndexes = find(diff(frameSignal)==1);

incrementSample = diff(frameIndexes);

incrementFrame = round(incrementSample/samplePerFrame);

errors =  find(incrementFrame>1);

for i = flip(1:1:numel(errors))

    j = errors(i);
    k = incrementFrame(j);
    insertion = 1:1:k-1; 
    insertion = insertion * round(incrementSample(j)/incrementFrame(j)) +frameIndexes(j);
    frameIndexes = [frameIndexes(1:j), insertion, frameIndexes(j+1:end)];

end

if numel(errors) > 0
    msgbox([num2str(numel(errors)) ' missing frames have been corrected.'])
    figure
    histogram(incrementFrame)
    title('nFrame/frameIndex')
    
end
    
end