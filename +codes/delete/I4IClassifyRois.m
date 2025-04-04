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



%   sData = vr.iman.classifyROIs(sData, sDataDir);

%save(fullfile(sDataDir,sessionID),'sData');
% close all
disp([sessionID ' - ROI classification done'])

save(fullfile(sDataDir,sessionID),'sData');
    

disp([ '- DONE - sData fle has been saved.'])

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