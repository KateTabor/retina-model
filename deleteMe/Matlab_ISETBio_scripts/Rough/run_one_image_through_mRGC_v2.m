% run_one_image_through_mRGC.m  (script version)

%% --- Init ISETBio/ISETCam 
isInit = false;
if exist('ieSessionGet','file') == 2
    try isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try evalin('base','ieInit;'); catch, ieInit; end
end

%% --- Fixed locations (exact paths you provided)
opticsFile = '/Users/kate/Documents/retina-model/optics/optics_humanWVF_3.0mm.mat';
lmFile     = '/Users/kate/Documents/retina-model/mosaics/mRGC_2deg_LM.mat';
lmsFile    = '/Users/kate/Documents/retina-model/mosaics/mRGC_2deg_LMS.mat';

%% --- Single image to process  *** CHANGED ***
thisImg = '/Users/kate/Documents/retina-model/bars_by1_b+100_XYZ.tiff';

% Save responses in the SAME directory as the image  *** CHANGED ***
[imgDir, imgName, ~] = fileparts(thisImg);
outDir = imgDir;   % no dated subfolder

assert(exist(thisImg,'file')==2, 'Image not found: %s', thisImg);

%% --- Load saved optics and mosaics
S = load(opticsFile, 'oi');            oi   = S.oi;
S = load(lmFile,   'theMRGCmosaic');   mLM  = S.theMRGCmosaic;
S = load(lmsFile,  'theMRGCmosaic_S'); mLMS = S.theMRGCmosaic_S;

%% --- Deterministic time settings
dt        = 0.010;                      % 10 ms
nFrames   = 10;                         % 10 frames
timeAxis  = (0:nFrames-1) * dt;

%% --- Process JUST THIS ONE IMAGE  *** CHANGED ***
outFiles = cell(1,1);
N  = 1;
fprintf('Processing %d image...\n', N);
t0 = tic;

ii = 1;
% (you could put this in a for-loop for future, but not needed for N=1)

% Build scene (no resize), modest FOV
sceneFOVdeg = 2;
scene = sceneFromXYZfloat32(thisImg, sceneFOVdeg);

% Optical image via pre-saved optics
oiNow = oiCompute(oi, scene);

% Run both mosaics (deterministic, avg over time)
respLM  = runOneMosaic(mLM,  oiNow, dt, nFrames, timeAxis);
respLMS = runOneMosaic(mLMS, oiNow, dt, nFrames, timeAxis);

% ----- Save simple output in SAME DIR as image  *** CHANGED ***
parent = imgDir;   % not really used in tag, but keep structure minimal
[~, parentName] = fileparts(parent);    % e.g., 'retina-model'
tag = sprintf('%s_%s_dt%gms_T%d', parentName, imgName, dt*1000, nFrames);

outFile       = fullfile(outDir, ['mRGCresp_' tag '.mat']);
rgc_positions = mLM.rgcRFpositionsDegs;

save(outFile, 'thisImg', 'respLM', 'respLMS', 'rgc_positions', ...
              'dt', 'nFrames', 'opticsFile', 'lmFile', 'lmsFile', '-v7');
outFiles{ii} = outFile;

elapsed = toc(t0);
fprintf('Done. Elapsed %.1fs. Saved to:\n  %s\n', elapsed, outFile);

%% ------------------------------------------------------------------------

