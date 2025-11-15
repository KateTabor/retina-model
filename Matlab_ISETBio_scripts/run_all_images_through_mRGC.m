function outFiles = run_all_images_through_mRGC()
% Process ALL images under .../image-set/{train,val,test}/trial_XXXXXX/.
% Deterministic: dt=10 ms, T=10 frames, average over time.
%
% No per-image stats or figures; keeps the progress meter.

%% --- Init ISETBio/ISETCam (your standard snippet)
isInit = false;
if exist('ieSessionGet','file') == 2
    try, isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try, evalin('base','ieInit;'); catch, ieInit; end
end

%% --- Fixed locations (exact paths you provided)
opticsFile = '/Users/kate/Documents/MATLAB/Toolboxes/isetbio/optics_humanWVF_3.0mm.mat';
lmFile     = '/Users/kate/Documents/MATLAB/Toolboxes/isetbio/mRGC_2deg_LM.mat';
lmsFile    = '/Users/kate/Documents/MATLAB/Toolboxes/isetbio/mRGC_2deg_LMS.mat';

%% --- Output root -> dated subfolder with abbreviated function name
outRoot = '/Users/kate/Documents/retina-model';
stamp   = datestr(now,'yyyymmdd');
outDir  = fullfile(outRoot, [stamp '_runAllImages']);   % e.g., 20251108_runAllImages
if ~exist(outDir,'dir'), mkdir(outDir); end

%% --- Load saved optics and mosaics
S = load(opticsFile, 'oi');            oi   = S.oi;
S = load(lmFile,   'theMRGCmosaic');   mLM  = S.theMRGCmosaic;
S = load(lmsFile,  'theMRGCmosaic_S'); mLMS = S.theMRGCmosaic_S;

%% --- Gather ALL XYZ.tiff images from train/val/test
rootDir = '/Users/kate/Documents/retina-model/image-set'; 
splits  = {'train','val','test'}; 
imgList = {};
for s = 1:numel(splits)
    dd = dir(fullfile(rootDir, splits{s}, '**', '*_XYZ.tif*'));
    for k = 1:numel(dd)
        imgList{end+1} = fullfile(dd(k).folder, dd(k).name); %#ok<AGROW>
    end
end
assert(~isempty(imgList), 'No *_XYZ.tif* files found under %s/{train,val,test}', rootDir);

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

    % Run both mosaics (deterministic, avg over time)
    respLM  = runOneMosaic(mLM,  oiNow, dt, nFrames, timeAxis);
    respLMS = runOneMosaic(mLMS, oiNow, dt, nFrames, timeAxis);

    % ----- Save simple output (to split subfolder; no stats/figs)
    [imgDir,imgName,~] = fileparts(thisImg);
    [~, parent]        = fileparts(imgDir);    % e.g., trial_000123
    tag                = sprintf('%s_%s_dt%gms_T%d', parent, imgName, dt*1000, nFrames);

    % Detect split from path (train/val/test) — no 'trial' here
    if contains(thisImg, [filesep 'train' filesep])
        splitName = 'train';
    elseif contains(thisImg, [filesep 'val' filesep])
        splitName = 'val';
    elseif contains(thisImg, [filesep 'test' filesep])
        splitName = 'test';
    else
        splitName = 'misc';
    end

    % Ensure split subfolder exists under the dated root
    splitOutDir = fullfile(outDir, splitName);
    if ~exist(splitOutDir,'dir'), mkdir(splitOutDir); end

    outFile       = fullfile(splitOutDir, ['mRGCresp_' tag '.mat']);
    rgc_positions = mLM.rgcRFpositionsDegs;

    save(outFile, 'thisImg', 'respLM', 'respLMS', 'rgc_positions', ...
                  'dt', 'nFrames', 'opticsFile', 'lmFile', 'lmsFile', '-v7');
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

%% ------------------------------------------------------------------------
function respMean = runOneMosaic(mosaic, oi, dt, nFrames, timeAxis)
% Deterministic mRGC run: use mosaic.inputConeMosaic, frozen noise, avg over time.

cm = mosaic.inputConeMosaic;     % confirmed on your install

% Cones: deterministic integration
cm.noiseFlag = 'frozen';
if isprop(cm, 'integrationTime'), cm.integrationTime = dt; end

% Also fix mRGC noise policy to deterministic (frozen noise)
if isprop(mosaic, 'noiseFlag'), mosaic.noiseFlag = 'frozen'; end

% Single-frame cone response to the OI
coneResp = cm.compute(oi);
coneResp = double(squeeze(coneResp));

% Reshape to [1 x T x Ncones] and replicate to nFrames
if isvector(coneResp)
    nCones = numel(coneResp);
    cone3D = reshape(coneResp, [1 1 nCones]);
else
    sz = size(coneResp);
    if numel(sz)==2 && all(sz > 1)
        if sz(1) > sz(2), tmp = coneResp.'; else, tmp = coneResp; end
        cone3D = reshape(tmp, [1 size(tmp,1) size(tmp,2)]);
    else
        if numel(sz) < 3, coneResp(:,:,2) = coneResp; sz = size(coneResp); end
        tmp    = reshape(coneResp, [], sz(3));
        cone3D = reshape(tmp.', [1 sz(3) numel(tmp)/sz(3)]);
    end
end
if size(cone3D,2) == 1
    cone3D = repmat(cone3D, [1 nFrames 1]);
else
    if size(cone3D,2) < nFrames
        cone3D = repmat(cone3D, [1 ceil(nFrames/size(cone3D,2)) 1]);
    end
    cone3D = cone3D(:,1:nFrames,:);
end

% Run mRGC mosaic (deterministic seed; do not pass nTrials/timeResolutionSeconds)
[rgcResp, rgcRespNoisy] = mosaic.compute(cone3D, timeAxis, 'seed', 1);

% Prefer noisy instances if present (with frozen noise it's deterministic)
if ~isempty(rgcRespNoisy)
    X = squeeze(rgcRespNoisy);   % [T x N]
else
    X = squeeze(rgcResp);        % [T x N]
end
if isvector(X), X = X(:)'; end
respMean = mean(X,1);            % [1 x N] average over time
end

%% ------------------------------------------------------------------------
function scene = sceneFromXYZfloat32(tiffPath, sceneFOVdeg)
% Minimal XYZ(D65) -> scene helper (linear sRGB, spectral display).

xyz = double(imread(tiffPath));            % HxWx3, float32->double
if size(xyz,3) ~= 3
    error('Expected HxWx3 XYZ image. Got size: %s', mat2str(size(xyz)));
end
if max(xyz(:)) > 2
    warning('XYZ seems >2.0; confirm writer scaling (expected ~0..1).');
end

% XYZ -> *linear* sRGB (no gamma)
if exist('colorTransformMatrix','file') == 2
    M_xyz2srgb = colorTransformMatrix('xyz2srgb');
else
    M_xyz2srgb = [ 3.2406 -1.5372 -0.4986; ...
                  -0.9689  1.8758  0.0415; ...
                   0.0557 -0.2040  1.0570 ];
end
sz = size(xyz);
rgb_lin = imageLinearTransform(reshape(xyz,[],3), M_xyz2srgb);
rgb_lin = reshape(rgb_lin, sz);
rgb_lin = ieClip(rgb_lin, 0, 1);

% Build spectral scene using display primaries (no gamma/uint8)
d     = displayCreate('LCD-Apple');
scene = sceneFromFile(rgb_lin, 'rgb', [], d);

% FOV
if nargin < 2 || isempty(sceneFOVdeg), sceneFOVdeg = 2; end
scene = sceneSet(scene, 'fov', sceneFOVdeg);
end
