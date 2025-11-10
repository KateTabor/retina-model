function theMRGCmosaic = build2degONmRGCMosaic_v2(varargin)
% build2degONmRGCMosaic  Load + visualize the 2-deg temporal ON mRGC mosaic (small figs).
%
% Syntax:
%   theMRGCmosaic = build2degONmRGCMosaic_v2
%   theMRGCmosaic = build2degONmRGCMosaic_v2('exemplarRGCindex', 10)
%
% Description:
%   (i)   loads the prebaked ON midget RGC mosaic centered ~2 deg temporal,
%   (ii)  visualizes the input cone mosaic (L/M/S cones),
%   (iii) visualizes the mRGC RF centers for the whole patch,
%   (iv)  visualizes center+surround cone pooling map for one exemplar RGC.
%
% Name-Value pairs:
%   'coneMosaicSpecies'         - 'human' (default) or 'macaque'
%   'opticsSubjectName'         - default: 'PLOSpaperDefaultSubject'
%   'rgcMosaicName'             - default: 'PLOSpaperTemporal2DegsMosaic'
%   'targetVisualSTFdescriptor' - default: 'default'
%   'exemplarRGCindex'          - index for detailed map (default: auto)
%   'closeOpenFigures'          - true/false (default: true)
%
% Returns:
%   theMRGCmosaic  - the loaded mRGC mosaic (for use in the workspace)

% --- keep caller args safe
userArgs = varargin;

% --- init ISETBio/ISETCam (avoid varargin being cleared)
isInit = false;
if exist('ieSessionGet','file') == 2
    try isInit = ieSessionGet('initialized'); catch, isInit = false; end
end
if ~isInit && exist('ieInit','file') == 2
    try evalin('base','ieInit;'); catch, ieInit; end
end

% --- parse inputs (simple)
p = inputParser;
p.addParameter('coneMosaicSpecies','human');
p.addParameter('opticsSubjectName','PLOSpaperDefaultSubject');
p.addParameter('rgcMosaicName','PLOSpaperTemporal2DegsMosaic');
p.addParameter('targetVisualSTFdescriptor','default');
p.addParameter('exemplarRGCindex',NaN);
p.addParameter('closeOpenFigures',true);
p.parse(userArgs{:});
opts = p.Results;

if opts.closeOpenFigures, close all; end

% --- load the prebaked mosaic + PSF (compute PSF for overlay if needed later)
[theMRGCmosaic, ~, ~] = mRGCMosaic.loadPrebakedMosaic( ...
    opts.coneMosaicSpecies, ...
    opts.opticsSubjectName, ...
    opts.rgcMosaicName, ...
    opts.targetVisualSTFdescriptor, ...
    'computeTheMosaicOptics', true);

% --- basic footprint printout
allPos = theMRGCmosaic.rgcRFpositionsDegs;  % [N x 2] deg
nRGCs  = theMRGCmosaic.rgcsNum;
minX = min(allPos(:,1));  maxX = max(allPos(:,1));
minY = min(allPos(:,2));  maxY = max(allPos(:,2));
fprintf('Loaded: %s | species=%s | subject=%s\n', ...
    opts.rgcMosaicName, opts.coneMosaicSpecies, opts.opticsSubjectName);
fprintf('Center [deg]=[%.2f %.2f] | Size [deg]=[%.2f %.2f] | RGCs=%d\n', ...
    theMRGCmosaic.eccentricityDegs(1), theMRGCmosaic.eccentricityDegs(2), ...
    theMRGCmosaic.sizeDegs(1), theMRGCmosaic.sizeDegs(2), nRGCs);
fprintf('X=[%.2f %.2f]  Y=[%.2f %.2f]\n', minX, maxX, minY, maxY);

% --- domain limits/ticks (centered on mosaic center)
visualizedWidthDegs  = theMRGCmosaic.sizeDegs(1);
visualizedHeightDegs = theMRGCmosaic.sizeDegs(2);
domainVisualizationLimits = [ ...
    theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs  * [-1 1], ...
    theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 1]];
domainVisualizationTicks = struct( ...
    'x', theMRGCmosaic.eccentricityDegs(1) + 0.5 * visualizedWidthDegs  * [-1 -0.5 0 0.5 1], ...
    'y', theMRGCmosaic.eccentricityDegs(2) + 0.5 * visualizedHeightDegs * [-1 -0.5 0 0.5 1]);

%% (1) CONE MOSAIC (L/M/S types)
cm = theMRGCmosaic.inputConeMosaic;   % confirmed on your install
hFig1 = figure(1); clf; set(hFig1,'Position',[200 200 640 420]); % small fig
cm.visualize( ...
    'figureHandle', hFig1, ...
    'domainVisualizationLimits', domainVisualizationLimits, ...
    'domainVisualizationTicks',  domainVisualizationTicks, ...
    'plotTitle', 'Input CONE mosaic (L/M/S types)');

%% (2) mRGC RF centers (center footprints)
minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap; % "Chichilnisky-style"
hFig2 = figure(2); clf; set(hFig2,'Position',[200 680 640 420]); % small fig
theMRGCmosaic.visualize( ...
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

% Use a simple, robust center-relative surround threshold (matches ISETBio examples)
minSurroundConeWeight = 1e-3;  % relative to the center

theMRGCmosaic.visualizeCenterSurroundConePoolingMap( ...
    targetIdx, ...
    'minConeWeightForVisualizingRFcenterPooling',   minCenterConeWeight, ...
    'minConeWeightForVisualizingRFsurroundPooling', minSurroundConeWeight, ...
    'minSurroundConeWeightRelativity', 'center', ...  % your build accepts {'center','surround'}
    'withLineWeightingFunctions', true, ...
    'scaleBarDegs', 0.1, ...
    'doNotLabelScaleBar', true, ...
    'plotTitle', sprintf('RGC #%d  (center + surround)', targetIdx), ...
    'figNo', 3, ...
    'figPos', [860 200], ...
    'withCustomFigureFormat', '1x1 giant rectangular-wide mosaic');

try set(figure(3), 'Position', [860 200 640 420]); end

save('mRGC_2deg_LM.mat','theMRGCmosaic','-v7.3');

end
