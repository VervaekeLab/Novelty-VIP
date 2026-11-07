%% STANDALONE RAW PUPIL TRACKER (MOUSE-FOLDER ISOLATED)
function pupilTrackerRaw()
    %% 1. SELECT RAW FILE & PARSE METADATA
    [rawFile, rawPath] = uigetfile('*.raw', 'Select the .raw Pupil Video');
    if isequal(rawFile, 0), return; end
    vidPath = fullfile(rawPath, rawFile);
    
    [vidWidth, vidHeight] = parseMetadata(rawPath);
    
    if isempty(vidWidth) || isempty(vidHeight)
        prompt = {'Enter Video Width (Pixels):', 'Enter Video Height (Pixels):'};
        dlgtitle = 'Metadata Not Found';
        dims = [1 50]; definput = {'250', '250'};
        answer = inputdlg(prompt, dlgtitle, dims, definput);
        if isempty(answer), return; end
        vidWidth = str2double(answer{1});
        vidHeight = str2double(answer{2});
    end
    
    fInfo = dir(vidPath);
    bytesPerFrame = vidWidth * vidHeight; 
    nFrames = floor(fInfo.bytes / bytesPerFrame);
    
    fprintf('Loaded raw file: %s\n', rawFile);
    fprintf('Resolution: %dx%d | Total Frames: %d\n', vidWidth, vidHeight, nFrames);

    %% 2. INTERNAL STATE
    state = struct('vidPath', vidPath, 'vidDir', rawPath, 'W', vidWidth, 'H', vidHeight, 'nFrames', nFrames, ...
                   'currentFrame', 1, 'threshold', 150, 'sensitivity', 0.95, ...
                   'roiCenter', [vidWidth/2, vidHeight/2], 'roiRadius', min(vidWidth, vidHeight)/3, ...
                   'pDiam', [], 'pX', [], 'pY', []);
                   
    %% 3. GUI SETUP
    uiFig = uifigure('Name', ['Raw Pupil Tracker (Bright Pupil) - ' rawFile], 'Position', [100 100 1000 650]);
    uiFig.Color = 'w';
    
    glMain = uigridlayout(uiFig, [2, 2]); 
    glMain.RowHeight = {'1x', 150}; glMain.ColumnWidth = {'1x', '1x'};
    
    axVid = uiaxes(glMain, 'BackgroundColor', 'k'); axVid.Title.String = 'Raw Video Preview';
    axMask = uiaxes(glMain, 'BackgroundColor', 'k'); axMask.Title.String = 'Threshold & Tracking Mask';
    disableDefaultInteractivity(axVid); disableDefaultInteractivity(axMask);
    
    pnlCtrl = uipanel(glMain, 'Title', 'Tracking Controls', 'BackgroundColor', 'w');
    pnlCtrl.Layout.Row = 2; pnlCtrl.Layout.Column = 1;
    glC = uigridlayout(pnlCtrl, [3, 2]); glC.ColumnWidth = {120, '1x'};
    
    lblFr = uilabel(glC, 'Text', 'Frame: 1', 'FontWeight', 'bold');
    sldFr = uislider(glC, 'Limits', [1 nFrames], 'Value', 1, 'MajorTicks', [], 'MinorTicks', [], 'ValueChangedFcn', @(s,e) updateSettings('frame', e.Value));
    
    lblTh = uilabel(glC, 'Text', 'Threshold: 150', 'FontWeight', 'bold');
    sldTh = uislider(glC, 'Limits', [0 255], 'Value', 150, 'MajorTicks', [], 'MinorTicks', [], 'ValueChangedFcn', @(s,e) updateSettings('thresh', e.Value));
    
    lblSe = uilabel(glC, 'Text', 'Sensitivity: 0.95', 'FontWeight', 'bold');
    sldSe = uislider(glC, 'Limits', [0.8 1.0], 'Value', 0.95, 'MajorTicks', [], 'MinorTicks', [], 'ValueChangedFcn', @(s,e) updateSettings('sens', e.Value));
    
    pnlAct = uipanel(glMain, 'Title', 'Processing', 'BackgroundColor', 'w');
    pnlAct.Layout.Row = 2; pnlAct.Layout.Column = 2;
    btnProcess = uibutton(pnlAct, 'Position', [20 50 180 40], 'Text', '1. Process & Compare', 'BackgroundColor', [0 .5 .2], 'FontColor', 'w', 'FontWeight', 'bold', 'ButtonPushedFcn', @processVideo);
    btnSave = uibutton(pnlAct, 'Position', [220 50 180 40], 'Text', '2. Save pupilData.mat', 'BackgroundColor', [0 0.45 0.74], 'FontColor', 'w', 'FontWeight', 'bold', 'Enable', 'off', 'ButtonPushedFcn', @saveData);

    roiObj = [];
    initPreview();

    %% 4. CALLBACKS
    function initPreview()
        fid = fopen(state.vidPath, 'r'); 
        fseek(fid, (state.currentFrame - 1) * bytesPerFrame, 'bof');
        img = fread(fid, [state.W, state.H], 'uint8')'; fclose(fid);
        
        imagesc(axVid, img, [0 255]); colormap(axVid, gray); axis(axVid, 'image'); axis(axVid, 'off');
        roiObj = drawcircle(axVid, 'Center', state.roiCenter, 'Radius', state.roiRadius, 'Color', 'y');
        addlistener(roiObj, 'MovingROI', @(src, ev) runTrackingPreview());
        runTrackingPreview();
    end

    function updateSettings(type, val)
        if strcmp(type, 'frame')
            state.currentFrame = round(val);
            lblFr.Text = sprintf('Frame: %d', state.currentFrame);
            initPreview();
        elseif strcmp(type, 'thresh')
            state.threshold = val;
            lblTh.Text = sprintf('Threshold: %d', round(val));
            runTrackingPreview();
        elseif strcmp(type, 'sens')
            state.sensitivity = val;
            lblSe.Text = sprintf('Sensitivity: %.2f', val);
            runTrackingPreview();
        end
    end

    function runTrackingPreview()
        if isempty(roiObj) || ~isvalid(axVid), return; end
        state.roiCenter = roiObj.Center; state.roiRadius = roiObj.Radius;
        
        fid = fopen(state.vidPath, 'r'); 
        fseek(fid, (state.currentFrame - 1) * bytesPerFrame, 'bof');
        img = fread(fid, [state.W, state.H], 'uint8')'; fclose(fid);
        
        [X, Y] = meshgrid(1:state.W, 1:state.H);
        mask = (X - state.roiCenter(1)).^2 + (Y - state.roiCenter(2)).^2 <= state.roiRadius^2;
        roiArea = img; 
        roiArea(~mask) = 0; 
        
        bin = imfill((roiArea > state.threshold), 'holes'); 
        imagesc(axMask, bin, [0 1]); colormap(axMask, gray); axis(axMask, 'image'); axis(axMask, 'off');
        
        hold(axVid, 'on'); delete(findobj(axVid, 'Tag', 'previewOverlay'));
        [c, r, m] = imfindcircles(bin, [10 round(state.roiRadius)], 'Sensitivity', state.sensitivity);
        if ~isempty(r)
            [~,i] = max(m); 
            hCirc = viscircles(axVid, c(i,:), r(i), 'Color', 'r', 'LineStyle', '--'); set(hCirc, 'Tag', 'previewOverlay');
            plot(axVid, c(i,1), c(i,2), 'r+', 'MarkerSize', 15, 'LineWidth', 2, 'Tag', 'previewOverlay');
        end
        drawnow limitrate;
    end

    %% 5. PROCESS ENTIRE VIDEO
    function processVideo(~, ~)
        d = uiprogressdlg(uiFig, 'Title', 'Tracking Pupil', 'Message', 'Processing video frames...', 'Cancelable', 'on');
        
        pDiam = NaN(state.nFrames, 1); pX = NaN(state.nFrames, 1); pY = NaN(state.nFrames, 1);
        
        c_roi = state.roiCenter; r_roi = state.roiRadius;
        x1 = max(1, round(c_roi(1)-r_roi)); x2 = min(state.W, round(c_roi(1)+r_roi));
        y1 = max(1, round(c_roi(2)-r_roi)); y2 = min(state.H, round(c_roi(2)+r_roi));
        [X, Y] = meshgrid(1:(x2-x1+1), 1:(y2-y1+1));
        circM = (X - (c_roi(1)-x1+1)).^2 + (Y - (c_roi(2)-y1+1)).^2 <= r_roi^2;
        
        fid = fopen(state.vidPath, 'r');
        for f = 1:state.nFrames
            if d.CancelRequested, break; end
            
            data = fread(fid, state.W * state.H, 'uint8');
            if isempty(data), break; end
            img = reshape(data, [state.W, state.H])'; 
            
            roi = img(y1:y2, x1:x2);
            roi(~circM) = 0; 
            
            bin = imfill((roi > state.threshold), 'holes'); 
            bin = bwareaopen(bin, 20);
            
            s = regionprops(bin, 'Area', 'EquivDiameter', 'Centroid', 'Eccentricity');
            if ~isempty(s)
                [~, id] = max([s.Area]); 
                if s(id).Eccentricity < 0.9 
                    pDiam(f) = s(id).EquivDiameter; 
                    pX(f) = s(id).Centroid(1) + x1 - 1; 
                    pY(f) = s(id).Centroid(2) + y1 - 1; 
                end
            end
            
            if mod(f, 1000) == 0, d.Value = f / state.nFrames; end
        end
        fclose(fid);
        
        state.pDiam = pDiam; state.pX = pX; state.pY = pY;
        btnSave.Enable = 'on';
        close(d);
        
        runDiagnosticComparison();
    end

    %% 6. DIAGNOSTIC COMPARISON
    function runDiagnosticComparison()
        matFiles = dir(fullfile(state.vidDir, '*.mat'));
        oldData = [];
        loadedFileName = '';
        
        for i = 1:length(matFiles)
            if startsWith(matFiles(i).name, '._'), continue; end
            
            try
                tmp = load(fullfile(matFiles(i).folder, matFiles(i).name));
                if isfield(tmp, 'Center') && isfield(tmp, 'Radius')
                    oldData = tmp;
                    loadedFileName = matFiles(i).name;
                    break;
                end
            catch
                continue;
            end
        end
        
        if isempty(oldData)
            uialert(uiFig, ['No baseline tracker data containing fields ''Center'' and ''Radius'' was found in: ' state.vidDir], 'Comparison Aborted');
            return;
        end
        
        fprintf('Comparing against baseline file: %s\n', loadedFileName);
        
        oldR = oldData.Radius(:);
        oldX = oldData.Center(:, 1);
        oldY = oldData.Center(:, 2);
        
        trimLen = min(state.nFrames, length(oldR));
        t = 1:trimLen;
        
        newD = state.pDiam(1:trimLen);
        newX = state.pX(1:trimLen);
        newY = state.pY(1:trimLen);
        
        oldX = oldX(1:trimLen);
        oldY = oldY(1:trimLen);
        
        newNaNs = sum(isnan(newD));
        oldNaNs = sum(isnan(oldR * 2));
        
        compFig = figure('Name', ['Diagnostics: ' loadedFileName], 'Color', 'w', 'Position', [150 150 1200 900]);
        tl = tiledlayout(3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
        
        titleStr = sprintf('Mouse Data Comparison | Reference File: %s\nTotal Frames: %d | Dropouts (NaNs) -> New Tracker: %d (%.1f%%) vs Baseline: %d (%.1f%%)', ...
            loadedFileName, trimLen, newNaNs, (newNaNs/trimLen)*100, oldNaNs, (oldNaNs/trimLen)*100);
        sgtitle(titleStr, 'FontWeight', 'bold', 'FontSize', 12, 'Interpreter', 'none');
        
        % Panel 1: Center Comparison (Fixed Color Strategy)
        ax1 = nexttile; hold(ax1, 'on');
        plot(ax1, t, oldX, 'Color', [0.65 0.85 1.0], 'LineWidth', 1.5, 'DisplayName', 'Baseline X (Light Blue)');
        plot(ax1, t, newX, 'Color', [0.00 0.45 0.74], 'LineWidth', 1.0, 'DisplayName', 'New X (Dark Blue)');
        plot(ax1, t, oldY, 'Color', [1.00 0.70 0.70], 'LineWidth', 1.5, 'DisplayName', 'Baseline Y (Light Red)');
        plot(ax1, t, newY, 'Color', [0.85 0.33 0.10], 'LineWidth', 1.0, 'DisplayName', 'New Y (Dark Red)');
        ylabel(ax1, 'Center Position (Pixels)'); legend(ax1, 'Location', 'best', 'NumColumns', 2);
        title(ax1, 'Pupil Coordinates Tracking');
        
        % Panel 2: Radius Comparison
        ax2 = nexttile; hold(ax2, 'on');
        plot(ax2, t, oldR, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.5, 'DisplayName', 'Baseline Radius');
        plot(ax2, t, newD/2, 'Color', [0.46 0.67 0.18], 'LineWidth', 1.0, 'DisplayName', 'New Radius');
        ylabel(ax2, 'Radius (Pixels)'); legend(ax2, 'Location', 'best');
        title(ax2, 'Pupil Size Resolution');
        
        % Panel 3: Stability/Derivatives Comparison
        ax3 = nexttile; hold(ax3, 'on');
        oldJitter = abs([0; diff(oldR)]);
        newJitter = abs([0; diff(newD/2)]);
        
        plot(ax3, t, oldJitter, 'Color', [0.75 0.75 0.75], 'LineWidth', 1.0, 'DisplayName', 'Baseline abs(\Delta Radius)');
        plot(ax3, t, newJitter, 'Color', [0.49 0.18 0.55], 'LineWidth', 1.0, 'DisplayName', 'New abs(\Delta Radius)');
        ylabel(ax3, 'Frame Delta (Pixels)'); xlabel(ax3, 'Frame Indices'); 
        legend(ax3, 'Location', 'best'); title(ax3, 'Tracking Stability (Derivative Spikes = Jumps/Errors)');
        
        linkaxes([ax1, ax2, ax3], 'x');
        xlim(ax1, [1 trimLen]);
    end

    %% 7. SAVE
    function saveData(~, ~)
        [file, path] = uiputfile('*.mat', 'Save Pupil Data', 'pupilData.mat');
        if isequal(file,0) || isequal(path,0), return; end 
        
        pupilData = struct('pupilDiam', state.pDiam, 'pupilX', state.pX, 'pupilY', state.pY);
        save(fullfile(path, file), 'pupilData');
        uialert(uiFig, 'Saved successfully.', 'Saved');
    end

    %% HELPER: PARSE METADATA
    function [w, h] = parseMetadata(folderPath)
        w = []; h = [];
        files = [dir(fullfile(folderPath, '*.ini')); dir(fullfile(folderPath, '*.yaml')); dir(fullfile(folderPath, '*.txt'))];
        for i = 1:length(files)
            try
                txt = fileread(fullfile(files(i).folder, files(i).name));
                tokW = regexp(txt, '(?i)(?:width|xpixels|x_size|x)[\s:=]+(\d+)', 'tokens');
                tokH = regexp(txt, '(?i)(?:height|ypixels|y_size|y)[\s:=]+(\d+)', 'tokens');
                if ~isempty(tokW) && isempty(w), w = str2double(tokW{1}{1}); end
                if ~isempty(tokH) && isempty(h), h = str2double(tokH{1}{1}); end
                if ~isempty(w) && ~isempty(h), return; end 
            catch
                continue;
            end
        end
    end
end