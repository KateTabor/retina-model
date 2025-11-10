function outFile = run_one_image_through_mRGC(imagePath)
% Minimal test driver: one image -> optics -> both mRGC mosaics -> save responses.
%
% Uses deterministic settings: dt=10 ms, T=10 frames, average over time.
%
% Output .mat is saved in /Users/kate/Documents/retina-model with a simple name.

%% --- Init ISETBio/ISETCam (your standard snippet)
isInit = false;
if exist('ieSessionGet','file') == 2
    try isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try evalin('base','ieInit;'); catch, ieInit; end
end

%% --- Fixed locations (exact paths you provided)
opticsFile = '/Users/kate/Documents/MATLAB/Toolboxes/isetbio/optics_humanWVF_3.0mm.mat';
lmFile     = '/Users/kate/Documents/MATLAB/Toolboxes/isetbio/mRGC_2deg_LM.mat';
lmsFile    = '/Users/kate/Documents/MATLAB/Toolboxes/isetbio/mRGC_2deg_LMS.mat';

%% --- Find a default image if none passed
if nargin < 1 || isempty(imagePath)
    rootDir = '/Users/kate/Documents/retina-model/image-set/train';
    dd = dir(fullfile(rootDir, '**', '*.tif*'));  % float32 XYZ (D65) TIFFs
    assert(~isempty(dd), 'No TIFFs found under %s', rootDir);
    imagePath = fullfile(dd(1).folder, dd(1).name);
end
assert(exist(imagePath,'file')==2, 'Image not found: %s', imagePath);

%% --- Load saved optics and mosaics
S = load(opticsFile, 'oi');               oi = S.oi;
S = load(lmFile,   'theMRGCmosaic');      mLM  = S.theMRGCmosaic;
S = load(lmsFile,  'theMRGCmosaic_S');    mLMS = S.theMRGCmosaic_S;

%% --- Build scene from float32 XYZ (D65) TIFF (no resizing), set modest FOV
sceneFOVdeg = 2;            % 2° matches mosaic framing
scene = sceneFromXYZfloat32(imagePath, sceneFOVdeg);

%% --- Compute optical image using the pre-saved optics
oi = oiCompute(oi, scene);

%% --- Deterministic time settings
dt        = 0.010;          % 10 ms
nFrames   = 10;             % 10 frames total
timeAxis  = (0:nFrames-1) * dt;

%% --- Run both mosaics and average over time
respLM  = runOneMosaic(mLM,  oi, dt, nFrames, timeAxis);
respLMS = runOneMosaic(mLMS, oi, dt, nFrames, timeAxis);

%% --- Save minimal output (simple/flat)
[imgDir,imgName,~] = fileparts(imagePath);
[~, parent]        = fileparts(imgDir);               % e.g., trial folder
tag                = sprintf('%s_%s_dt%gms_T%d', parent, imgName, dt*1000, nFrames);
outDir             = '/Users/kate/Documents/retina-model';
if ~exist(outDir,'dir'), mkdir(outDir); end
outFile            = fullfile(outDir, ['mRGCresp_' tag '.mat']);

% Include positions once (from LM), bookkeeping, and your paths
rgc_positions = mLM.rgcRFpositionsDegs;
save(outFile, 'imagePath', 'respLM', 'respLMS', 'rgc_positions', ...
              'dt', 'nFrames', 'opticsFile', 'lmFile', 'lmsFile', '-v7');
fprintf('Saved: %s\n', outFile);
end

%% ------------------------------------------------------------------------
function respMean = runOneMosaic(mosaic, oi, dt, nFrames, timeAxis)
% Deterministic mRGC run: get cone mosaic, set noise frozen, build 10 frames, avg over time.

% 1) find internal cone mosaic handle (robust to naming)
cm = [];
for nm = {'inputConeMosaic','coneMosaic','theConeMosaic','theInputConeMosaic'}
    if isprop(mosaic, nm{1}), cm = mosaic.(nm{1}); break; end
end
assert(~isempty(cm), 'No internal cone mosaic found in mRGC mosaic.');

% 2) deterministic settings on cones
cm.noiseFlag = 'frozen';           % deterministic noise
if isprop(cm, 'integrationTime'), cm.integrationTime = dt; end

% 3) single-frame cone response to the OI
coneResp = cm.compute(oi);         % returns cone samples (implementation-dependent)
coneResp = double(squeeze(coneResp));

% 4) reshape to [1 x T x Ncones] and replicate to nFrames
if isvector(coneResp)
    nCones = numel(coneResp);
    cone3D = reshape(coneResp, [1 1 nCones]);
else
    sz = size(coneResp);
    if numel(sz)==2 && all(sz > 1)
        % [nCones x T] or [T x nCones]; prefer [T x nCones]
        if sz(1) > sz(2), tmp = coneResp.'; else, tmp = coneResp; end
        cone3D = reshape(tmp, [1 size(tmp,1) size(tmp,2)]);
    else
        % [R x C x T] -> [1 x T x Ncones]
        if numel(sz) < 3, coneResp(:,:,2) = coneResp; sz = size(coneResp); end
        tmp    = reshape(coneResp, [], sz(3));
        cone3D = reshape(tmp.', [1 sz(3) numel(tmp)/sz(3)]);
    end
end
if size(cone3D,2) == 1
    cone3D = repmat(cone3D, [1 nFrames 1]);
else
    % If more than one frame came back, pad/trim to nFrames
    if size(cone3D,2) < nFrames
        cone3D = repmat(cone3D, [1 ceil(nFrames/size(cone3D,2)) 1]);
    end
    cone3D = cone3D(:,1:nFrames,:);
end

% 5) run mRGC mosaic (prefer noisy response if available)
[rgcResp, rgcRespNoisy] = mosaic.compute(cone3D, timeAxis, 'seed', 1);

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
% Minimal XYZ(D65) -> scene helper for our pipeline.
% Reads 512x512 float32 XYZ, converts to *linear* sRGB (no gamma/uint8),
% and builds a spectral scene via a calibrated display model.
% (Assumes ieInit already called by the caller.)

% 1) read float32 XYZ
xyz = double(imread(tiffPath));            % HxWx3, float32 -> double
if size(xyz,3) ~= 3
    error('Expected HxWx3 XYZ image. Got size: %s', mat2str(size(xyz)));
end
if max(xyz(:)) > 2
    warning('XYZ seems >2.0; confirm writer scaling (expected ~0..1).');
end

% 2) XYZ -> *linear* sRGB  (no gamma)
if exist('colorTransformMatrix','file') == 2
    M_xyz2srgb = colorTransformMatrix('xyz2srgb');
else
    % IEC 61966-2-1 linear RGB from XYZ (D65)
    M_xyz2srgb = [ 3.2406 -1.5372 -0.4986; ...
                  -0.9689  1.8758  0.0415; ...
                   0.0557 -0.2040  1.0570 ];
end
sz = size(xyz);
rgb_lin = imageLinearTransform(reshape(xyz,[],3), M_xyz2srgb);
rgb_lin = reshape(rgb_lin, sz);
rgb_lin = ieClip(rgb_lin, 0, 1);           % keep linear space

% 3) Build spectral scene using display primaries (no gamma, no uint8)
d     = displayCreate('LCD-Apple');        % spectral display model
scene = sceneFromFile(rgb_lin, 'rgb', [], d);

% 4) set FOV
if nargin < 2 || isempty(sceneFOVdeg), sceneFOVdeg = 2; end
scene = sceneSet(scene, 'fov', sceneFOVdeg);
end
