
%  
% 


function sData = classifyROIs(sData,sDataDir)


nFOVs =  length(sData.imdata);

for fov = 1:1:nFOVs
    
    sData = classifyROIsFOV(sData, fov, sDataDir);
end

end


function sData = classifyROIsFOV(sData, fov, sDataDir)


    SMI = []; SMIcorr = [];
    for t = 1:1:length(sData.trials.contextsMeta)
        trials = sData.trials.contextsMeta(t).trials;
        dataM = sData.imdata.binnedRoisDff(trials,:,:);

        [~, ~, SMIcorr(:,t)] = codes.getSMICorrAllRois(dataM(1:2:numel(trials),:,:),dataM(2:2:numel(trials),:,:),0.95);

        nROIs = sData.imdata.nROIs;

        for r = 1:1:nROIs
            SMI(r,t) = codes.getSMI(dataM(:,:,r));
        end

    end

    for r = 1:1:nROIs

        sData.imdata.roiMeta(r).smiPeak = SMI(r,:);
        sData.imdata.roiMeta(r).smiCorr = SMIcorr(r,:);

    end



%% Find place fields

smoothWindow = 5;
avThr = 0.05;
%singThr = 0.15;
stdThrPeak = 2.5;
stdThrOnset = 2;
reliabilityThr = 0.6;
movWindow = 5; 
fieldWidthThr = [3, 40]; % i.e. 6 - 80 cm
maxFieldNumber = 3;

nROIs = sData.imdata.nROIs;


for r = 1:1:nROIs
    
    nFields = [];
    
    for b = 1:1:length(sData.trials.contextsMeta)
        % Generate dataM from the trial blocks and +- trials to test
        % reliability
        trials = sData.trials.contextsMeta(b).trials;
        lastBlock = false;
        if max(trials) + movWindow*2-1 <= sData.trials.contextsMeta(end).blockEnd && min(trials) - (movWindow-1) > 0
            relTrials = min(trials)- (movWindow-1):1:(max(trials) + movWindow*2-1); relTrials = relTrials - (trials(1)-1);           
            trials = min(trials)- (movWindow-1):1:(max(trials) + movWindow*2-1);
        elseif b == 1 % first block
            relTrials = min(trials):1:(max(trials) + movWindow*2-1);           
            trials = min(trials):1:(max(trials) + movWindow*2-1);
        elseif max(trials) + movWindow*2-1 > sData.trials.contextsMeta(end).blockEnd % last block
            relTrials = min(trials)- (movWindow-1):1:max(trials); relTrials = relTrials - (trials(1)-1);           
            trials = min(trials)- (movWindow-1):1:max(trials);
            lastBlock = true;
        end
        
%         if max(trials) + movWindow-1 <= sData.trials.contextsMeta(end).blockEnd
%             trials = min(trials):1:(max(trials) + movWindow-1);
%         else
%             trials = min(trials):1:max(trials);
%         end
        
        dataM = sData.imdata.binnedRoisDff(trials,:,r);
        
        if sData.trials.contextsMeta(b).contextIndicator > 10
            ctxInd = [sData.trials.contextsMeta(:).contextIndicator];
            st = min([sData.trials.contextsMeta(ctxInd>0).blockStart]); en = max([sData.trials.contextsMeta(ctxInd>0).blockEnd]);
            dataCtx = sData.imdata.binnedRoisDff(st:en,:,r);
            dataM = dataCtx;
        else
            dataCtx = dataM;
        end
        
        % % Shift the largest peak to the center of the track
        [~, maxInd] = max(nanmean(dataM));
        sz = size(dataM); szc = size(dataCtx);
        dataM = reshape(circshift(reshape(dataM',1,[]),(65-maxInd)), [sz(2), sz(1)])';
        dataCtx = reshape(circshift(reshape(dataCtx',1,[]),(65-maxInd)), [szc(2), szc(1)])';
        trialIndMatr = ones(size(dataM)); for i = 1:1:size(trialIndMatr,1); trialIndMatr(i,:) = trialIndMatr(i,:)*i; end
        trialIndMatr = reshape(circshift(reshape(trialIndMatr',1,[]),(65-maxInd)), [sz(2), sz(1)])';
        if maxInd < 65
            trialIndMatr(1,1:(65-maxInd)) = NaN;
            dataM(1,1:(65-maxInd)) = NaN;
        end
        shiftXax = circshift(1:1:125,(65-maxInd));
        
        % Norm to max act of roi
        peakActivity = max(max(sData.imdata.binnedRoisDff(:,:,r))); % quantile(max(sData.imdata.binnedRoisDff(:,:,r),[],2),0.95) sData.imdata.roiMeta(r).peakDff;
        %minActivity =         quantile(sData.imdata.binnedRoisDff(:,:,r),0.05,2)
        dataM = smoothdata(dataM,2,'gaussian',smoothWindow)/peakActivity; %sData.imdata.roiMeta(r).peakDff;
        dataM = dataM - quantile(dataM,0.05,2);
        dataCtx = smoothdata(dataCtx,2,'gaussian',smoothWindow)/peakActivity; %sData.imdata.roiMeta(r).peakDff;
        dataCtx = dataCtx - quantile(dataCtx,0.05,2);
        
        
        
        
        % peakActThrVsStd(r,:) = [(peakActivity)*singThr, std(dataM(:),"omitnan")];
        
        %prominThr = (peakActivity)*singThr;
        %prominThr = std(dataM(:),"omitnan")*3;
        prominThr = std(reshape(smoothdata(sData.imdata.binnedRoisDff(:,:,r),2,'gaussian',smoothWindow)/peakActivity,1,[]),"omitnan")*stdThrPeak;

 %%%%%%%%       
        if ~lastBlock
            normAvData = nanmean(dataCtx(end-19-(movWindow*2-1):end-(movWindow*2-1),:)); 
        else
            normAvData = nanmean(dataCtx(end-19:end,:)); 
        end
        
        % Find potential peaks
        [~,avLocs, width, ~] = findpeaks([normAvData(end-49:end), normAvData, normAvData(1:50)],'MinPeakProminence',avThr);
        
        avLocs = avLocs - 50; % handle circularity
        ignore = avLocs < 1 | avLocs > 125;
        avLocs = avLocs(~ignore); width = width(~ignore);
        
        fieldData = struct([]);
        for field = 1:1:numel(avLocs)
           % isActive = false(1,numel(trials));
            peakLocInTr = nan(1,numel(trials));
            onsetLocInTr = nan(1,numel(trials));
            peakAmplInTr = nan(1,numel(trials));
            
            fieldRange = [floor(avLocs(field)-width(field)), ceil(avLocs(field)+width(field))];
            indNotDetect = true; endNotDetect = true;
            indTrial = NaN;
            fieldData(field).indTrial = NaN; fieldData(field).indPos = NaN; fieldData(field).indSpeed = NaN; 
            fieldData(field).fieldEndTrial = NaN; fieldData(field).reliability = NaN; fieldData(field).fieldWidth = NaN;
            for t = 1:1:numel(trials)
                % Correct for circularity
                if t > 1 && t < size(dataM,1)
                    trialData = [dataM(t-1,end-49:end), dataM(t,:), dataM(t+1,1:50)];
                    %[~,singLocs] = findpeaks([dataM(t-1,end-49:end), dataM(t,:), dataM(t+1,1:50)],'MinPeakProminence',prominThr); % max(pks)*singThr %(max(dataM(:))-min(dataM(:)))      
                elseif t == 1
                    trialData = [dataM(t,end-49:end), dataM(t,:), dataM(t+1,1:50)];
                    %[~,singLocs] = findpeaks([dataM(t,end-49:end), dataM(t,:), dataM(t+1,1:50)],'MinPeakProminence',prominThr); % max(pks)*singThr
                else
                    trialData = [dataM(t-1,end-49:end), dataM(t,:), dataM(t,1:50)];
                    %[~,singLocs] = findpeaks([dataM(t-1,end-49:end), dataM(t,:), dataM(t,1:50)],'MinPeakProminence',prominThr); % max(pks)*singThr
                end
                [~,singLocs] = findpeaks(trialData,'MinPeakProminence',prominThr); % max(pks)*singThr
                
                onset = []; % for BTSP PF onset
                if numel(singLocs) > 0
                    for l = 1:1:numel(singLocs)
                        tempOnset = max(find((trialData(1:singLocs(l)) < prominThr/stdThrPeak*stdThrOnset)));
                        if numel(tempOnset) > 0
                            onset = [onset, max(find((trialData(1:singLocs(l)) < prominThr/stdThrPeak*stdThrOnset)))];
                        else
                            onset = [onset, singLocs(l)];
                        end
                    end
                end
                
                singLocs = singLocs - 50;
                onset = onset - 50;
                % ignore = singLocs < 1 | singLocs > 125; singLocs = singLocs(~ignore);
        
%                 if any(fieldRange> 125 | fieldRange < 1)
%                      singLocs((singLocs-avLocs(field)) < 1)
%                 end
                inField = fieldRange(1)< singLocs &  singLocs < fieldRange(2);
                singLocs(singLocs > 125) = 125; singLocs(singLocs < 1) = 1;
                onset(onset > 125) = 125; onset(onset < 1) = 1;
                tNonShift = min(trialIndMatr(t,singLocs(inField)));
                if isnan(tNonShift); inField(:) = false; tNonShift = t; end % only happens at the very beginning og a block
                % isActive(tNonShift) = any(inField);  %use the smallest trial index if there are multiple peaks
                if ~any(inField)
                    peakLocInTr(tNonShift) = NaN;
                    onsetLocInTr(tNonShift) = NaN;
                    peakAmplInTr(tNonShift) = NaN;
                else
                    peakLocInTr(tNonShift) = shiftXax(floor(mean(singLocs(inField))));                   
                    onsetLocInTr(tNonShift) = shiftXax(floor(mean(onset(inField))));
                    tempData = smoothdata(sData.imdata.binnedRoisDff(trials(tNonShift),:,r),2,'gaussian',smoothWindow);
                    peakAmplInTr(tNonShift) =  tempData(peakLocInTr(tNonShift));
                end
            end
            
            % Find induction trial and end trial. Here everything is real, non shifted trials and positions! 
            
            for t = 1:1:length(peakLocInTr)-(movWindow-1)
                if ~isnan(peakLocInTr(t)) && mean(~isnan(peakLocInTr(t:t+movWindow-1))) >= reliabilityThr && indNotDetect && t <= (length(peakLocInTr)-(movWindow*2-1))
                    indTrial = relTrials(t); %trialIndMatr(t,floor(mean(singLocs(inField))));
                        indTrialIndex = t; % this is the index within dataM
                    fieldData(field).indTrial = indTrial;
                    fieldData(field).indPos = peakLocInTr(t);
                    fieldData(field).indPosOnset = onsetLocInTr(t);
                    vel = smoothdata(sData.behavior.trialMatrices.binVel(trials(t),:),2,'gaussian',5);
                    fieldData(field).fieldEndTrial = length(peakLocInTr); %will be updated if an earlier field end is detected
                    
                    if peakLocInTr(t) < 5
                        fieldData(field).indSpeed = mean(vel(1:peakLocInTr(t))); % single bin!
                    else
                        fieldData(field).indSpeed = mean(vel(peakLocInTr(t)-4:peakLocInTr(t))); % mean of previous 5 bin
                    end
                    indNotDetect = false;
                elseif indNotDetect
                    peakLocInTr(t) = NaN;
                    onsetLocInTr(t) = NaN;
                end
                if (mean(~isnan(peakLocInTr(t:t+movWindow-1))) < reliabilityThr && ~indNotDetect && endNotDetect) ...
                        || (((length(peakLocInTr)-(movWindow-1))==t) && endNotDetect && ~indNotDetect ) % end of the field if reliability drops below threshold, or last trial (subrtact mov window from the beginning)
                    
                    fieldEndTrial = relTrials(t) + movWindow-2; % one trial back
                    endTrialIndex = t + movWindow-2;
                    peakLocInTr(endTrialIndex+1:end) = NaN;
                        onsetLocInTr(endTrialIndex+1:end) = NaN;
                        peakAmplInTr(endTrialIndex+1:end) = NaN;
                    notNan = find(~isnan(peakLocInTr));
%                     try
                    fieldData(field).fieldEndTrial = notNan(end)-(movWindow-1); % fieldEndTrial; % subtract the movWindow from the beginning!
                    fieldData(field).fieldDuration = fieldData(field).fieldEndTrial - fieldData(field).indTrial;
                    %                     catch
%                         fieldData(field).fieldEndTrial = fieldEndTrial; % ;
%                     end
                    endNotDetect = false;
                end

            end

            
            % ?
            % if indTrial > numel(trials)-5 || isnan(indTrial);  indTrial = numel(trials)-5; end % consider at least 5 trials for reliability!
            if ~isnan(indTrial)
            avData = mean(smoothdata(sData.imdata.binnedRoisDff(trials(indTrialIndex:endTrialIndex),:,r),2,'gaussian',smoothWindow));
            fieldData(field).peakPos = shiftXax(avLocs(field));
            fieldData(field).peakAmpl = avData(shiftXax(avLocs(field))) - quantile(avData,0.05); % from induction
            fieldData(field).fieldWidth = width(field); % from induction
            
%             peakLocInTr = peakLocInTr(movWindow:end);
%             	peakAmplInTr = peakAmplInTr(movWindow:end);
            peakLocInTr(indTrialIndex) = fieldData(field).indPos;
            peakLocInTr = peakLocInTr - fieldData(field).indPos;
            peakLocInTr(peakLocInTr > 62.5) = -62.5 + (peakLocInTr(peakLocInTr > 62.5) - 62.5);
            peakLocInTr(peakLocInTr <-62.5) = 62.5 + (peakLocInTr(peakLocInTr <-62.5) + 62.5);
            
            onsetLocInTr(indTrialIndex) = fieldData(field).indPosOnset;
            onsetLocInTr = onsetLocInTr - fieldData(field).indPosOnset;
            onsetLocInTr(onsetLocInTr > 62.5) = -62.5 + (onsetLocInTr(onsetLocInTr > 62.5) - 62.5);
            onsetLocInTr(onsetLocInTr <-62.5) = 62.5 + (onsetLocInTr(onsetLocInTr <-62.5) + 62.5);

            if b > 1
                fieldData(field).peakLocInTr = peakLocInTr(movWindow:end);
                fieldData(field).onsetLocInTr = onsetLocInTr(movWindow:end);
                fieldData(field).peakAmplInTr = peakAmplInTr(movWindow:end);
            else
                fieldData(field).peakLocInTr = peakLocInTr(1:end);
                fieldData(field).onsetLocInTr = onsetLocInTr(1:end);
                fieldData(field).peakAmplInTr = peakAmplInTr(1:end);
            end

            fieldData(field).reliability = mean(~isnan(peakLocInTr(indTrialIndex:endTrialIndex)));% from induction
            %if ~isfield(fieldData(field),'indPos') || fieldData(field).indTrial > numel(trials)-5
            %   fieldData(field).peakShift = NaN;
            %else
            fieldData(field).peakShift = nanmean(peakLocInTr(indTrialIndex+1:endTrialIndex));% from induction
            fieldData(field).onsetShift = nanmean(onsetLocInTr(indTrialIndex+1:endTrialIndex));% from induction
            %end
            % fieldData(field).isSignificant = fieldData(field).reliability >= reliabilityThr;
            end          
            
        end
        
        
        lowCorr = sData.imdata.roiMeta(r).smiPeak(b) < 2 | sData.imdata.roiMeta(r).smiCorr(b) < 2;
        lowAct = false; % sData.imdata.roiMeta(r).activityLevel < 0.05; sum([sData.imdata.roiMeta(:).activityLevel] < 0.05);
        lowSNR = sData.imdata.roiMeta(r).SNR < 2;
        
        % Delete non-significant fields
        if ~isempty( fieldData)
            keepFields = [fieldData(:).reliability] >= reliabilityThr & ...
                [fieldData(:).fieldWidth] > fieldWidthThr(1) & ...
                [fieldData(:).fieldWidth] < fieldWidthThr(2) & ...
               [fieldData(:).fieldEndTrial] > 1 & ...
...                [fieldData(:).fieldEndTrial] >= movWindow & ...
                ([fieldData(:).fieldEndTrial] - [fieldData(:).indTrial]) >= movWindow & ...
                true;
            if lowAct || lowCorr || lowSNR || sum(keepFields) > maxFieldNumber
                keepFields(:) = false;
            end
            fieldData = fieldData(keepFields);
        end
             
        
        sData.imdata.roiMeta(r).placeFieldsInBlock(b).fieldData = fieldData;
        nFields(b) = length(fieldData);
        
    end
    
%     for b = 2:1:4
%     
%     
%     
%     end
    
    %nFields(lowCorr) = 0; if lowAct; nFields = nFields*0; end
    
    sData.imdata.roiMeta(r).nFields = nFields;
    sData.imdata.roiMeta(r).placeCell = nFields > 0;
    
end


% update imdata.placeCells

placeCells = reshape([sData.imdata.roiMeta(:).placeCell], length(sData.trials.contextsMeta),[]);

blockData = struct();

for b = 1:1:length(sData.trials.contextsMeta)
    
    blockData(b).placeCells = find(placeCells(b,:));
    
    blockData(b).nPlaceCells = numel(find(placeCells(b,:)));
    blockData(b).placeCellFraction = numel(find(placeCells(b,:)))/nROIs;
    
    
end

sData.imdata.placeCells = blockData;




%% Control Figure

%{

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



%}




end

