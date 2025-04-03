For analyisis Matlab R2018a or R2021a was used.
Installing Matlab 2021a: download the installer for your platform and MATLAB release from MathWorks Downloads, installation takes 30-45 minutes.
System requirements: 
https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2021a-windows.pdf
https://www.mathworks.com/content/dam/mathworks/mathworks-dot-com/support/sysreq/files/system-requirements-release-2021a-macintosh.pdf

*************************************************************************************************************************************
What does the sData file contain?

Information about the session:
sData.sessionInfo

Information about the animal, surgery details:
sData.mouseInfo

Raw behavioral data:
sData.daqdata: recorded by LabView
The sampling rate is sData.daqdata.meta.fs = 3000. 

Trial blocks of contexts and optogenetic-stimulation protocol: 
sData.trials.trialBlocksMeta: contains information about the trial block structure

Behavioral and opto-stimulation data:
sData.behavior: contain information about the animal position on the wheel, animal speed, if the animal licks, when the reward was given, when the optogenetic stimulation was applied.

Imaging data: 
Calcium signals were recorded with SciScan and ROIs were detected and signal was extracted using NANSEN (https://github.com/VervaekeLab/NANSEN, Hennestad et al., submitted).
sData.imdata.roiSignals
sData.imdata.roiSignals(2).roif: raw Ca signals for each RIOs (rows) and each frame (columns)
sData.imdata.roiSignals(2).npilf: raw neurophil signals around each RIOs (rows) and each frame (columns)
sData.imdata.roiSignals(2).dff: dF/F signals for each RIOs (rows) and each frame (columns) calculated by NANSEN (neurophil is subtracted)
sData.imdata.roiSignals(2).deconv: deconvolved dF/F signals for each RIOs (rows) and each frame (columns) calculated by NANSEN (using CalmAn)
(For all analysis the sData.imdata.roiSignals(2).dff was used.) 

sData.imdata also contain information about recording parameters (sData.imdata.meta), signal extraction options (sData.imdata.signalExtractionOptions),
ROI statistics (sData.imdata.roiStat).

************************************************************************************************************************************
TEST SESSIONS

/example sessions 
Folder contains example sessions for a VIP-ArchT, VIP-ChrimsonR and control experiment. 
************************************************************************************************************************************
ANALYSIS PIPELINE (all codes were written by Mate Neubrandt)

Run 'pipeline.mat' to analyze imaging data in the selected session (sData file).
Expected run time for a standard desktop computer: each session takes 12-16 minutes.

The resulting data will be added as extra fields of the sData structures such as:

sData.imdata.roiStat: basic statistics (SNR, activity level)
sData.imdata.roiMeta: contains for each roi the basic statistics, result of place cell/ place field classification
sData.imdata.binnedRoisDff: 3D structure of the binned imaging data: 1) trials, 2) track position, 3) ROI index
sData.imdata.avBinnedRois: position tuning curves for each ROIs, i.e. average activity of all trials within each trial block 
sData.imdata.placeCells: place cell data within each trial block

sData.behavior.signals: I downsampled the behavior data to the imaging sampling rate, which is sData.imdata.meta.fps = ~31. 