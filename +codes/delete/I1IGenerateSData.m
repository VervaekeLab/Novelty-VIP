function varargout = GenerateSData(sessionObject, varargin)
%GENERATESDATA Summary of this function goes here
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
    ATTRIBUTES = {'batch', 'unqueueable'};   

    
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
    


% select single or multiple files to analyze 
[fileNames,filePath,~] = uigetfile('*.tdms','Select ALL!!! TDMS files from a LabView recording folder','C:\','MultiSelect','on' );


if ~iscell(fileNames)
    fileNames = {fileNames};
end
    
fileNumber = size(fileNames,2);




sDataDir = strcmp({sessionObject(1).DataLocation.Name}, 'sData'); 
sDataDir = sessionObject(1).DataLocation(sDataDir).RootPath;

tempFolderRoot = 'C:\TempData';
dateFolder = filePath(end-10:end);
tempFolder = fullfile(tempFolderRoot,dateFolder);

if ~isfolder(tempFolder); mkdir(tempFolder); end


fileStruct = dir(filePath);
%targetPaths(1:length(fileStruct)) = {''};
for f = 1:1:length(fileStruct)
    if ~fileStruct(f).isdir % exclude '.' & '..'

        sourceFile = fullfile(fileStruct(f).folder,fileStruct(f).name);
        destinationFile = fullfile(tempFolder,fileStruct(f).name);
        
        copyfile(sourceFile,destinationFile)

    end
end



%% Run analysis and generate sData files


tic;
for f = 1:1:fileNumber

    fileName = fileNames{f};
    try
        vr.sData.AnalyseRawData_V11_Remap(fileName,tempFolder);
        disp([num2str(f) ' / ' num2str(fileNumber) ': ' fileName ' - Done'])
    catch
        disp([num2str(f) ' / ' num2str(fileNumber) ': ' fileName ' - ERROR (probably too short or interrupted session)'])
    end

end
toc;


%% Copy files to target folders

fileStruct = dir(tempFolder);


for f = 3:1:length(fileStruct)
    if numel(fileStruct(f).name) > 2 % exclude '.' & '..'

        mouseName = fileStruct(f).name(1:5);
        if params.NansenDefNaming
            targetPath = [sDataDir '\subject-' mouseName];
        else
            targetPath = [sDataDir '\' mouseName];
        end

        if ~isfolder(targetPath); mkdir(targetPath); end
        
        sourceFile = fullfile(fileStruct(f).folder,fileStruct(f).name);
        destinationFile = fullfile(targetPath,fileStruct(f).name);
        
        [~,~,ext] = fileparts(sourceFile);
        if ~strcmp(ext,'.tdms')
            copyfile(sourceFile,destinationFile)
        end

        delete(sourceFile)

    end
end
disp('Files have been copied to the mouse folder.')


    
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
    params.NansenDefNaming = false;

end