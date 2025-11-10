function theMRGCmosaic_S = build2degONmRGCMosaic_enableSconeSurrounds_v2(varargin)
% build2degONmRGCMosaic_enableSconeSurrounds_v2
% Enable S-cone input to the surrounds of ALL mRGCs and generate the same
% figures/plots as build2degONmRGCMosaic_v2, plus connectivity plots.
%
% Syntax:
%   theMRGCmosaic_S = build2degONmRGCMosaic_enableSconeSurrounds_v2
%   theMRGCmosaic_S = build2degONmRGCMosaic_enableSconeSurrounds_v2(spec{:}, 'exemplarRGCindex', 10)
%
% Name-Value pairs (mirrors v2 wherever possible):
%   'coneMosaicSpecies'         - 'human' (default) or 'macaque')
%   'opticsSubjectName'         - default: 'PLOSpaperDefaultSubject'
%   'rgcMosaicName'             - default: 'PLOSpaperTemporal2DegsMosaic'
%   'targetVisualSTFdescriptor' - default: 'default'
%   'exemplarRGCindex'          - index for detailed map (default: auto)
%   'indicesOfSconeSurroundPoolingEnabledRGCs' - [] (default = ALL RGCs)
%   'closeOpenFigures'          - true/false (default: true)
%
% Returns:
%   theMRGCmosaic_S  - modified mosaic with S-cone surrounds enabled

% --- keep caller args safe (consistent with v2)
userArgs = varargin;

% --- init ISETBio/ISETCam (same style as v2)
isInit = false;
if exist('ieSessionGet','file') == 2
    try isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try evalin('base','ieInit;'); catch, ieInit; end
end

% --- parse inputs (mirrors v2 + S-cone extras)
p = inputParser;
p.addParameter('coneMosaicSpecies','human');
p.addParameter('opticsSubjectName','PLOSpaperDefaultSubject');
p.addParameter('rgcMosaicName','PLOSpaperTemporal2DegsMosaic');
p.addParameter('targetVisualSTFdescriptor','default');
p.addParameter('exemplarRGCindex',NaN);
p.addParameter('indicesOfSconeSurroundPoolingEnabledRGCs',[]);
p.addParameter('closeOpenFigures',true);
p.addParameter('previewBeforeAfter', false);   % show before/after for 1–2 cells
p.addParameter('previewRGCindices', []);       % indices to preview (optional)
p.parse(userArgs{:});
opts = p.Results;

if opts.closeOpenFigures, close all; end

% --- load the SAME prebaked mosaic spec as v2 (optics not required here)
[theMRGCmosaic_S, ~, ~] = mRGCMosaic.loadPrebakedMosaic( ...
    opts.coneMosaicSpecies, ...
    opts.opticsSubjectName, ...
    opts.rgcMosaicName, ...
    opts.targetVisualSTFdescriptor, ...
    'computeTheMosaicOptics', false);

% --- basic footprint printout (same style as v2)
allPos = theMRGCmosaic_S.rgcRFpositionsDegs;  % [N x 2] deg
nRGCs  = theMRGCmosaic_S.rgcsNum;
minX = min(allPos(:,1));  maxX = max(allPos(:,1));
minY = min(allPos(:,2));  maxY = max(allPos(:,2));
fprintf('Loaded (S-enabled target): %s | species=%s | subject=%s\n', ...
    opts.rgcMosaicName, opts.coneMosaicSpecies, opts.opticsSubjectName);
fprintf('Center [deg]=[%.2f %.2f] | Size [deg]=[%.2f %.2f] | RGCs=%d\n', ...
    theMRGCmosaic_S.eccentricityDegs(1), theMRGCmosaic_S.eccentricityDegs(2), ...
    theMRGCmosaic_S.sizeDegs(1), theMRGCmosaic_S.sizeDegs(2), nRGCs);
fprintf('X=[%.2f %.2f]  Y=[%.2f %.2f]\n', minX, maxX, minY, maxY);

% --- domain limits/ticks (centered on mosaic center) - same code as v2
visualizedWidthDegs  = theMRGCmosaic_S.sizeDegs(1);
visualizedHeightDegs = theMRGCmosaic_S.sizeDegs(2);
domainVisualizationLimits = [ ...
    theMRGCmosaic_S.eccentricityDegs(1) + 0.5 * visualizedWidthDegs  * [-1 1], ...
    theMRGCmosaic_S.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 1]];
domainVisualizationTicks = struct( ...
    'x', theMRGCmosaic_S.eccentricityDegs(1) + 0.5 * visualizedWidthDegs  * [-1 -0.5 0 0.5 1], ...
    'y', theMRGCmosaic_S.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 -0.5 0 0.5 1]);

%% ADDED (optional before/after preview for 1–2 cells using a clean copy)
if opts.previewBeforeAfter
    if ~isempty(opts.previewRGCindices)
        previewIdx = opts.previewRGCindices(:)';
    else
        if isnan(opts.exemplarRGCindex)
            previewIdx = min(max(1, round(nRGCs/2)), nRGCs);
        else
            previewIdx = min(max(1, round(opts.exemplarRGCindex)), nRGCs);
        end
    end
    [mosaicPreview, ~, ~] = mRGCMosaic.loadPrebakedMosaic( ...
        opts.coneMosaicSpecies, opts.opticsSubjectName, opts.rgcMosaicName, ...
        opts.targetVisualSTFdescriptor, 'computeTheMosaicOptics', false);

    mosaicPreview.enableSconeSurroundPoolingInSelectCells( ...
        previewIdx, ...
        'visualizeBeforeAndAfterConePoolingMaps', true, ...
        'exportVisualizationPDFdirectory', '' );
end


% ===============================
% (A) ENABLE S-cone surrounds
% ===============================
% --- choose target RGCs (default = ALL)
if isempty(opts.indicesOfSconeSurroundPoolingEnabledRGCs)
    targetRGCs = 1:nRGCs;
else
    targetRGCs = opts.indicesOfSconeSurroundPoolingEnabledRGCs(:)';
end

fprintf('Enabling S-cone surround pooling for %d / %d mRGCs ...\n', ...
    numel(targetRGCs), nRGCs);

theMRGCmosaic_S.enableSconeSurroundPoolingInSelectCells( ...
    targetRGCs, ...
    'visualizeBeforeAndAfterConePoolingMaps', false, ...
    'exportVisualizationPDFdirectory', '' );
fprintf('Done. Centers unchanged; surrounds now pool indiscriminately from L/M/S.\n');

% ===============================
% (B) FIGURES — SAME AS v2
% ===============================

%% (1) CONE MOSAIC (L/M/S types)
cm = theMRGCmosaic_S.inputConeMosaic;
hFig1 = figure(1); clf; set(hFig1,'Position',[200 200 640 420]); % small fig
cm.visualize( ...
    'figureHandle', hFig1, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks',  domainVisualizationTicks, ...
    'plotTitle', 'Input CONE mosaic (L/M/S types)');

%% (2) mRGC RF centers (center footprints)
minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap; % "Chichilnisky-style"
hFig2 = figure(2); clf; set(hFig2,'Position',[200 680 640 420]); % small fig
theMRGCmosaic_S.visualize( ...
    'figureHandle', hFig2, ...
    'identifyInputCones', false, ...
    'identifyPooledCones', false, ...
    'identifiedConeAperture','lightCollectingArea4sigma', ...
    'identifiedConeApertureThetaSamples',16, ...
    'centerSubregionContourSamples',32, ...
    'minConeWeightVisualized', minCenterConeWeight, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks',  domainVisualizationTicks, ...
    'plotTitle', sprintf('mRGC RF centers (min center weight = %0.3f)', minCenterConeWeight));

%% (3) One exemplar: center + surround pooling maps
if isnan(opts.exemplarRGCindex)
    targetIdx = min(max(1, round(nRGCs/2)), nRGCs);
else
    targetIdx = min(max(1, round(opts.exemplarRGCindex)), nRGCs);
end

% Keep the same tutorial-style surround thresholding paradigm as v2
minSurroundConeWeight = 1e-3;   % relative to the center
theMRGCmosaic_S.visualizeCenterSurroundConePoolingMap( ...
    targetIdx, ...
    'minConeWeightForVisualizingRFcenterPooling',   minCenterConeWeight, ...
    'minConeWeightForVisualizingRFsurroundPooling', minSurroundConeWeight, ...
    'minSurroundConeWeightRelativity', 'center', ...   % IMPORTANT: allowed values {'center','surround'}
    'withLineWeightingFunctions', true, ...
    'scaleBarDegs', 0.1, ...
    'doNotLabelScaleBar', true, ...
    'plotTitle', sprintf('RGC #%d  (center + surround; S-enabled)', targetIdx), ...
    'figNo', 3, ...
    'figPos', [860 200], ...
    'withCustomFigureFormat', '1x1 giant rectangular-wide mosaic');
try, set(figure(3), 'Position', [860 200 640 420]); end

% ===============================
% (C) CONNECTIVITY PLOTS (like analyze_mRGCMosaic)
% ===============================

% Matrices are [nCones x nRGC] on your build
C = theMRGCmosaic_S.rgcRFcenterConeConnectivityMatrix;
S = theMRGCmosaic_S.rgcRFsurroundConeConnectivityMatrix;
pos = theMRGCmosaic_S.rgcRFpositionsDegs; % [nRGC x 2]

% --- Center threshold (figure-style)
centerThr = mRGCMosaic.sensitivityAtPointOfOverlap;

% --- Surround threshold: 0.001 x center peak, per RGC (tutorial style)
alpha = 1e-3;
centerPeaks = full(max(abs(C),[],1));       % [1 x nRGC], ensure full
thrPerCell  = alpha * centerPeaks;          % [1 x nRGC]

% --- Counts per RGC (sum across cones)
centerCounts   = full(sum(abs(C) >  centerThr, 1)).';                            % [nRGC x 1]
surroundCounts = full(sum(abs(S) >  bsxfun(@times, ones(size(S,1),1), thrPerCell), 1)).';  % [nRGC x 1]

% --- x-eccentricity
xEcc = pos(:,1);

% Plot: center cones per RGC vs x-ecc
figure(10); clf; set(gcf,'Position',[900 200 600 400]);
plot(xEcc, centerCounts, '.', 'MarkerSize', 6);
xlabel('x-eccentricity (deg, ~-2 \rightarrow 0)'); ylabel('# cones in CENTER');
title('CENTER: # pooled cones per mRGC vs x-ecc (S-enabled surrounds)');
grid on;

% Plot: surround cones per RGC vs x-ecc
figure(11); clf; set(gcf,'Position',[900 640 600 400]);
plot(xEcc, surroundCounts, '.', 'MarkerSize', 6);
xlabel('x-eccentricity (deg, ~-2 \rightarrow 0)'); ylabel('# cones in SURROUND');
title('SURROUND: # pooled cones per mRGC vs x-ecc (S-enabled surrounds)');
grid on;

% Quick text summary (avoid sparse fprintf issue by ensuring full)
fprintf('CENTER cones/RGC (S-enabled):   min=%d  median=%.1f  mean=%.2f  max=%d\n', ...
    min(centerCounts), median(centerCounts), mean(centerCounts), max(centerCounts));
fprintf('SURROUND cones/RGC (S-enabled): min=%d  median=%.1f  mean=%.2f  max=%d\n', ...
    min(surroundCounts), median(surroundCounts), mean(surroundCounts), max(surroundCounts));

save('mRGC_2deg_LMS.mat','theMRGCmosaic_S','-v7.3');

end
