function outFiles = run_all_images_through_mRGC_v3()
% Process ALL images under
% Deterministic: dt=10 ms, T=10 frames, average over time.
% keeps the progress meter.

%% --- Init ISETBio/ISETCam 
isInit = false;
if exist('ieSessionGet','file') == 2
    try isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try evalin('base','ieInit;'); catch, ieInit; end
end

%% --- Fixed locations (exact paths you provided)
opticsFile   = '/Users/kate/Documents/retina-model/optics/optics_humanWVF_3.0mm.mat';
lmFile       = '/Users/kate/Documents/retina-model/mosaics/mRGC_2deg_LM.mat';
lmsFile      = '/Users/kate/Documents/retina-model/mosaics/mRGC_2deg_LMS.mat';
normBitsFile = '/Users/kate/Documents/retina-model/mosaics/mRGC_2deg_LMS_normBits.mat';  %%% NEW

% --- Date tokens used by both input and output paths
today = datetime('today','Format','yyyy-MM-dd');   % e.g., 2025-11-13
stamp = char(datetime(today,'Format','yyyyMMdd')); % e.g., 20251113

% Input: the generator's dated folder 
rootDir = fullfile('/Users/kate/Documents/retina-model/2025-11-18_bars_gradients');

% Output: the response pass (your existing naming, same date)
outDir  = fullfile('/Users/kate/Documents/retina-model/2025-11-18_bars_gradients', ...
                   [stamp '_runAllImages']);
if ~exist(outDir,'dir'), mkdir(outDir); end

%% --- Load saved optics and mosaics
S  = load(opticsFile, 'oi');            oi   = S.oi;
S  = load(lmFile,   'theMRGCmosaic');   mLM  = S.theMRGCmosaic;
S  = load(lmsFile,  'theMRGCmosaic_S'); mLMS = S.theMRGCmosaic_S;
Sbits    = load(normBitsFile, 'LMSnorm');                                  %%% NEW
LMSnorm  = Sbits.LMSnorm;                                                  %%% NEW
S_override_norm = LMSnorm.S_surround;                                      %%% NEW

%% --- Gather ALL XYZ.tiff images 
imgList = {};
dd = dir(fullfile(rootDir, '*_XYZ.tif*'));
for k = 1:numel(dd)
    imgList{end+1} = fullfile(dd(k).folder, dd(k).name); 
end

assert(~isempty(imgList), 'No *_XYZ.tif* files found under %s', rootDir);

%% --- Deterministic time settings
dt        = 0.010;                      % 10 ms
nFrames   = 10;                         % 10 frames
timeAxis  = (0:nFrames-1) * dt;

%% --- Loop images (all) with progress
outFiles = cell(1,numel(imgList));
N  = numel(imgList);
fprintf('Processing %d images...\n', N);
t0 = tic;

for ii = 1:N
    thisImg = imgList{ii};
    assert(exist(thisImg,'file')==2, 'Image not found: %s', thisImg);

    % Build scene (no resize), modest FOV
    sceneFOVdeg = 2;
    scene = sceneFromXYZfloat32(thisImg, sceneFOVdeg);

    % Optical image via pre-saved optics
    oiNow = oiCompute(oi, scene);

    % Run LM and LMS mosaics (deterministic, avg over time)
    respLM      = runOneMosaic(mLM,  oiNow, dt, nFrames, timeAxis);
    respLMS     = runOneMosaic(mLMS, oiNow, dt, nFrames, timeAxis);
    respLMS_norm = runOneMosaic(mLMS, oiNow, dt, nFrames, timeAxis, S_override_norm);  %%% NEW

    % ----- Save simple output (to split subfolder; no stats/figs)
    [imgDir,imgName,~] = fileparts(thisImg);
    [~, parent]        = fileparts(imgDir);    % e.g., trial_000123
    tag                = sprintf('%s_%s_dt%gms_T%d', parent, imgName, dt*1000, nFrames);

    outFile       = fullfile(outDir, ['mRGCresp_' tag '.mat']);
    rgc_positions = mLM.rgcRFpositionsDegs;

    save(outFile, ...
         'thisImg', 'respLM', 'respLMS', 'respLMS_norm', ...          %%% CHANGED (added respLMS_norm)
         'rgc_positions', 'dt', 'nFrames', ...
         'opticsFile', 'lmFile', 'lmsFile', 'normBitsFile', '-v7');   %%% CHANGED (added normBitsFile)
    outFiles{ii} = outFile;

    % ---- Progress (every 10 images, plus first/last)
    if mod(ii,10)==0 || ii==1 || ii==N
        pct     = 100*ii/N;
        elapsed = toc(t0);
        rate    = ii/max(elapsed, eps);     % images/sec
        eta     = (N-ii)/max(rate, eps);    % seconds remaining
        fprintf('Progress: %d/%d (%.1f%%) | elapsed %.1fs | ETA %.1fs\r', ...
                ii, N, pct, elapsed, eta);
        drawnow limitrate nocallbacks
        if ii==N, fprintf('\n'); end
    end
end
end


