% analyze_mRGCMosaic.m
% Run after:
%   theMRGCmosaic = build2degONmRGCMosaic_v2;

%% 0) Sanity
if ~exist('theMRGCmosaic','var')
    error('Run first: theMRGCmosaic = build2degONmRGCMosaic_v2;');
end

%% 1) Cone counts in the INPUT cone mosaic
cm = theMRGCmosaic.inputConeMosaic;
nL = numel(cm.lConeIndices);
nM = numel(cm.mConeIndices);
nS = numel(cm.sConeIndices);
fprintf('Cone counts: L=%d, M=%d, S=%d, total=%d\n', nL, nM, nS, cm.conesNum);

%% 2) Center/surround pooling matrices and positions ([nCones x nRGC])
C   = theMRGCmosaic.rgcRFcenterConeConnectivityMatrix;   % [nCones x nRGC]
S   = theMRGCmosaic.rgcRFsurroundConeConnectivityMatrix; % [nCones x nRGC]
pos = theMRGCmosaic.rgcRFpositionsDegs;                  % [nRGC x 2], deg
nRGCs = theMRGCmosaic.rgcsNum;

if size(C,2) ~= nRGCs
    error('Expected C to be [nCones x nRGC]. Got %dx%d with nRGC=%d.', size(C,1), size(C,2), nRGCs);
end

% Thresholds (tutorial-style)
minCenterConeWeight = mRGCMosaic.sensitivityAtPointOfOverlap;  % scalar
alpha = 1e-3;                                                   % surround = 0.001 x center peak (per cell)

% Per-RGC center peak, then per-RGC surround threshold
centerPeaks = max(abs(C), [], 1);           % [1 x nRGC]
thrPerCell  = alpha * centerPeaks;          % [1 x nRGC]

% Counts per RGC (sum across cones)
centerCounts   = sum(abs(C) > minCenterConeWeight, 1).';  % [nRGC x 1], may be sparse
surroundCounts = sum(abs(S) > thrPerCell,          1).';  % [nRGC x 1], may be sparse

% >>> Make sure we use full vectors for stats/printing/plotting
centerCounts   = full(centerCounts);
surroundCounts = full(surroundCounts);

% x-eccentricity (temporal side negative ~ -2 → 0)
xEcc = pos(:,1);

% Quick stats
fprintf('CENTER cones/RGC:   min=%d  median=%.1f  mean=%.2f  max=%d\n', ...
    min(centerCounts), median(centerCounts), mean(centerCounts), max(centerCounts));
fprintf('SURROUND cones/RGC: min=%d  median=%.1f  mean=%.2f  max=%d (alpha=%g)\n', ...
    min(surroundCounts), median(surroundCounts), mean(surroundCounts), max(surroundCounts), alpha);

%% 3) Plot: # center cones per RGC vs x-ecc
figure(10); clf; set(gcf,'Position',[900 200 600 400]);
plot(xEcc, centerCounts, '.', 'MarkerSize', 6);
xlabel('x-eccentricity (deg, ~-2 \rightarrow 0)'); ylabel('# cones in CENTER');
title('CENTER: # pooled cones per mRGC vs x-ecc'); grid on;
xlim([min(xEcc) max(xEcc)]);

%% 4) Plot: # surround cones per RGC vs x-ecc
figure(11); clf; set(gcf,'Position',[900 640 600 400]);
plot(xEcc, surroundCounts, '.', 'MarkerSize', 6);
xlabel('x-eccentricity (deg, ~-2 \rightarrow 0)'); ylabel('# cones in SURROUND');
title(sprintf('SURROUND: # pooled cones per mRGC vs x-ecc (\\alpha=%g)', alpha)); grid on;
xlim([min(xEcc) max(xEcc)]);

%% 5) Optional histograms
figure(12); clf; set(gcf,'Position',[300 200 560 400]);
histogram(centerCounts, 100); grid on;
xlabel('# cones in CENTER'); ylabel('count of mRGCs'); title('CENTER counts distribution');

figure(13); clf; set(gcf,'Position',[300 640 560 400]);
histogram(surroundCounts, 100); grid on;
xlabel('# cones in SURROUND'); ylabel('count of mRGCs'); title(sprintf('SURROUND counts distribution (\\alpha=%g)', alpha));
