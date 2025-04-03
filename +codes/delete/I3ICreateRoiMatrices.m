function varargout = CreateRoiMatrices(sessionObject, varargin)
%CREATEROIMATRICES Summary of this function goes here
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
    ATTRIBUTES = {'serial', 'queueable'};   

    
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


for i = 1:1:length(sessionObject)    
try
fovLocation = {sessionObject(i).FovLocation};
catch
    fovLocation = {'N/A'};
end
if strcmp(fovLocation{1},'N/A')
% Manual input of fov location
prompt = {'FOV location 1','FOV location 2 (mesoscope)','FOV location 3 (mesoscope)','FOV location 4 (mesoscope)','FOV location 5 (mesoscope)','FOV location 6 (mesoscope)','FOV location 7 (mesoscope)','FOV location 8 (mesoscope)'};
dlgtitle = ['Specify Imaging Location for ' sessionObject(i).sessionID ' ! WARNING: if multiple session is selected this will apply to all of them!'];
definput = {'','','','','','','',''};
fovLocation = inputdlg(prompt,dlgtitle,[1 40],definput);
fovLocation = fovLocation(~strcmp(fovLocation,''));
end



ind = find(strcmp({sessionObject(i).DataLocation.Name},'sData')); mousePath = fullfile(sessionObject(i).DataLocation(ind).RootPath,sessionObject(i).subjectID);

try
    sDataDir = sessionObject(i).getSessionFolder('sData');
    file = vr.findFileInFolder(sDataDir,sessionObject(i).sessionID,'mat');
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


sData = vr.transposeRoiSignals(sData);


sData = vr.downsampleBehavior(sData); % save(fullfile(sDataDir,sessionID),'sData');
disp([sessionID ' - Behavior downsampled'])
%% quality check
try
if ~isfield(sData.imdata(1),'roiStat')
    sData = vr.iman.imQualityChecks(sData, sDataDir);
end
disp([sessionID ' - QC and ROI stats done'])
catch
disp('QC ERROR!')
end
if ~isfield(sData.imdata(1),'binnedRoisDff')
    sData = vr.iman.createRoiMatricesContexts(sData);
    save(fullfile(sDataDir,sessionID),'sData');
end
disp([sessionID ' - ROI matrices has been created'])


for fov = 1:1:length(sData.imdata)
    sData.imdata(fov).fovLocation = fovLocation{fov};
end



%{
sData = vr.iman.classifyROIs(sData, sDataDir);
%save(fullfile(sDataDir,sessionID),'sData');
% close all
disp([sessionID ' - ROI classification done'])
%}    
save(fullfile(sDataDir,sessionID),'sData');


disp([ '- DONE - sData file has been saved.'])

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