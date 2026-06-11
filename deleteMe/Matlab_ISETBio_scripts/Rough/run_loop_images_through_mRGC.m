function outFiles = run_loop_images_through_mRGC(imagePath)
% Minimal driver: image(s) -> optics -> both mRGC mosaics -> save responses.
% Deterministic: dt=10 ms, T=10 frames, average over time.

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

%% --- Output root -> dated subfolder with abbreviated function name
outRoot = '/Users/kate/Documents/retina-model';
stamp   = datestr(now,'yyyymmdd');
outDir  = fullfile(outRoot, [stamp '_runLoopImages']);   % e.g., 20251108_runLoopImages
if ~exist(outDir,'dir'), mkdir(outDir); end

%% --- Load saved optics and mosaics
S = load(opticsFile, 'oi');                   oi   = S.oi;
S = load(lmFile,   'theMRGCmosaic');          mLM  = S.theMRGCmosaic;
S = load(lmsFile,  'theMRGCmosaic_S');        mLMS = S.theMRGCmosaic_S;

%% --- Figure out images to process
imgList = {};
if nargin < 1 || isempty(imagePath)
    % take first 5 XYZ TIFFs in train split
    rootDir = '/Users/kate/Documents/retina-model/image-set/train';
    dd = dir(fullfile(rootDir, '**', '*.tif*'));
    assert(~isempty(dd), 'No TIFFs found under %s', rootDir);
    take = min(5, numel(dd));
    for k = 1:take, imgList{end+1} = fullfile(dd(k).folder, dd(k).name); end %#ok<AGROW>
else
    imgList = {imagePath};   % single test image
end

%% --- Deterministic time settings
dt        = 0.010;                      % 10 ms
nFrames   = 10;                         % 10 frames
timeAxis  = (0:nFrames-1) * dt;

%% --- Loop images (1 or 5)
outFiles = cell(1,numel(imgList));
% ---- One combined figure for all images
f = figure('Name','LM vs LMS: all images','NumberTitle','off');
ax = axes('Parent',f); hold(ax,'on');
cmap = lines(numel(imgList));  % distinct colors per image
for ii = 1:numel(imgList)
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

    % Save simple output
    [imgDir,imgName,~] = fileparts(thisImg);
    [~, parent]        = fileparts(imgDir);    % e.g., trial folder
    tag                = sprintf('%s_%s_dt%gms_T%d', parent, imgName, dt*1000, nFrames);
    outFile            = fullfile(outDir, ['mRGCresp_' tag '.mat']);
    rgc_positions      = mLM.rgcRFpositionsDegs;

    save(outFile, 'thisImg', 'respLM', 'respLMS', 'rgc_positions', ...
                  'dt', 'nFrames', 'opticsFile', 'lmFile', 'lmsFile', '-v7');
    fprintf('Saved: %s\n', outFile);
    outFiles{ii} = outFile;

    % -------- Stats (for now, print per image)
    numRGCs = size(rgc_positions,1);
    fprintf('Dims: LM=%d  LMS=%d  Ncells=%d\n', numel(respLM), numel(respLMS), numRGCs);
    fprintf('LM : min=%0.4g  max=%0.4g  mean=%0.4g\n',  min(respLM),  max(respLM),  mean(respLM));
    fprintf('LMS: min=%0.4g  max=%0.4g  mean=%0.4g\n\n', min(respLMS), max(respLMS), mean(respLMS));

    % ---- Plot this image onto the shared axes
    label = sprintf('%s/%s', parent, imgName);
    scatter(ax, respLM, respLMS, 6, cmap(ii,:), 'Marker','.', 'DisplayName', label);
end
xlabel(ax,'LM mean response'); ylabel(ax,'LMS mean response');
title(ax, sprintf('Per-cell mean responses (%d images)', numel(imgList)));
xlim(ax,[-10 25]); ylim(ax,[-10 25]); axis(ax,'square'); grid(ax,'on');
legend(ax,'Location','eastoutside');
drawnow;
end

%% ------------------------------------------------------------------------
function respMean = runOneMosaic(mosaic, oi, dt, nFrames, timeAxis)
% Deterministic mRGC run: use mosaic.theInputConeMosaic, frozen noise, avg over time.

% 1) we now **know** the property name
cm = mosaic.inputConeMosaic;

% Cones: deterministic integration
cm.noiseFlag = 'frozen';
if isprop(cm, 'integrationTime'), cm.integrationTime = dt; end

% Also fix mRGC noise policy to deterministic (frozen noise)
if isprop(mosaic, 'noiseFlag'), mosaic.noiseFlag = 'frozen'; end

% 2) single-frame cone response to the OI
coneResp = cm.compute(oi);         % size varies by cm config
coneResp = double(squeeze(coneResp));

% 3) reshape to [1 x T x Ncones] and replicate to nFrames
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

% 4) run mRGC mosaic (deterministic seed).  Do **not** pass 'nTrials' or 'timeResolutionSeconds'.
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
% Minimal XYZ(D65) -> scene helper for our pipeline (linear sRGB, spectral display).

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
