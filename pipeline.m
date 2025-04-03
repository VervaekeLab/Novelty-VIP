%% Pipeline

% Select and load sData file
[sessionID,sDataDir,~] = uigetfile('*.mat','Select sData File','','MultiSelect','off' );
load(fullfile(sDataDir,sessionID));
fileparts(fileparts(sDataDir))

savePath = fullfile(sDataDir,'analyzed'); 

%% downsample behavior data to match with imaging frames

sData = codes.downsampleBehavior(sData); % save(fullfile(sDataDir,sessionID),'sData');
disp([sessionID ' - Behavior downsampled'])
   
% save sData file
save(fullfile(savePath,['analyzed_' sessionID]),'sData');
disp([ '- DONE - sData file has been saved.'])

%% quality check

sData = codes.imQualityChecks(sData);
disp([sessionID ' - QC and ROI stats done'])

% save sData file
save(fullfile(savePath,['analyzed_' sessionID]),'sData');
disp([ '- DONE - sData file has been saved.'])

%% bin imaging data


sData = codes.createRoiMatricesContexts(sData);
disp([sessionID ' - ROI matrices has been created'])

% save sData file
save(fullfile(savePath,['analyzed_' sessionID]),'sData');
disp([ '- DONE - sData file has been saved.'])

%% classify ROIs

sData = codes.classifyROIs(sData,sDataDir);
disp([sessionID ' - place field have been detected'])

% save sData file
save(fullfile(savePath,['analyzed_' sessionID]),'sData');
disp([ '- DONE - sData file has been saved.'])
