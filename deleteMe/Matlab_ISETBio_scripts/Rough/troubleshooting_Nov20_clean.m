%% troubleshooting_Nov20_clean.m
% Assumes the following exist in the workspace:
%   mLM, mLMS                             : LM-only and LMS mosaic objects
%   S_surround_LMS_norm                   : normalized LMS surround connectivity [nCones_LMS x nRGC]
%   W_L_LMS_norm, W_M_LMS_norm, W_S_LMS_norm : normalized LMS L-, M-, and S-surround weights (optional)
%   stats                                 : struct from normalize_LMS_surround_to_LM_v2 (not used here)

if ~exist('mLM','var') || ~exist('mLMS','var')
    error('mLM and mLMS must be in the workspace.');
end
if ~exist('S_surround_LMS_norm','var')
    error('S_surround_LMS_norm not found. Run normalize_LMS_surround_to_LM_v2.m first.');
end

%% 0. Basic connectivity, cone mosaics, indices, positions

% Surround connectivity matrices
S_surround_LM  = mLM.rgcRFsurroundConeConnectivityMatrix;   % [nCones_LM  x nRGC]
S_surround_LMS = mLMS.rgcRFsurroundConeConnectivityMatrix;  % [nCones_LMS x nRGC]

% Cone mosaics
cm_LM  = mLM.inputConeMosaic;
cm_LMS = mLMS.inputConeMosaic;

% Cone-type indices (per mosaic)
lIdx_LM = cm_LM.lConeIndices;
mIdx_LM = cm_LM.mConeIndices;
sIdx_LM = cm_LM.sConeIndices;

lIdx_LMS = cm_LMS.lConeIndices;
mIdx_LMS = cm_LMS.mConeIndices;
sIdx_LMS = cm_LMS.sConeIndices;

% RGC positions
pos_LM  = mLM.rgcRFpositionsDegs;    % [nRGC x 2]
pos_LMS = mLMS.rgcRFpositionsDegs;

nRGC = size(pos_LM,1);
if size(pos_LMS,1) ~= nRGC
    error('LM and LMS mosaics must have the same number of RGCs.');
end

% Shared axis limits (LM + LMS mosaics)
xAll  = [pos_LM(:,1); pos_LMS(:,1)];
yAll  = [pos_LM(:,2); pos_LMS(:,2)];
xlims = [min(xAll) max(xAll)];
ylims = [min(yAll) max(yAll)];

%% 1. Cone-class surround weights per RGC (sum over cones)

% --- LM model weights ---
W_L_LM = full(S_surround_LM(lIdx_LM,:));      % [nL_LM x nRGC]
W_M_LM = full(S_surround_LM(mIdx_LM,:));      % [nM_LM x nRGC]
W_S_LM = full(S_surround_LM(sIdx_LM,:));      % [nS_LM x nRGC]

Lsum_LM = sum(W_L_LM, 1);                     % [1 x nRGC]
Msum_LM = sum(W_M_LM, 1);
Ssum_LM = sum(W_S_LM, 1);

% --- LMS model weights (original) ---
W_L_LMS = full(S_surround_LMS(lIdx_LMS,:));   % [nL_LMS x nRGC]
W_M_LMS = full(S_surround_LMS(mIdx_LMS,:));   % [nM_LMS x nRGC]
W_S_LMS = full(S_surround_LMS(sIdx_LMS,:));   % [nS_LMS x nRGC]

Lsum_LMS = sum(W_L_LMS, 1);
Msum_LMS = sum(W_M_LMS, 1);
Ssum_LMS = sum(W_S_LMS, 1);

% --- LMS model weights (normalized) ---
% Use precomputed W_L_LMS_norm / W_M_LMS_norm / W_S_LMS_norm if available;
% otherwise derive from S_surround_LMS_norm.
if exist('W_L_LMS_norm','var')
    W_L_LMS_norm = full(W_L_LMS_norm);
else
    W_L_LMS_norm = full(S_surround_LMS_norm(lIdx_LMS,:));
end

if exist('W_M_LMS_norm','var')
    W_M_LMS_norm = full(W_M_LMS_norm);
else
    W_M_LMS_norm = full(S_surround_LMS_norm(mIdx_LMS,:));
end

if exist('W_S_LMS_norm','var')
    W_S_LMS_norm = full(W_S_LMS_norm);
else
    W_S_LMS_norm = full(S_surround_LMS_norm(sIdx_LMS,:));
end

Lsum_LMS_norm = sum(W_L_LMS_norm, 1);
Msum_LMS_norm = sum(W_M_LMS_norm, 1);
Ssum_LMS_norm = sum(W_S_LMS_norm, 1);

% Common color scale for weight maps (include normalized LMS)
allW  = [Ssum_LM Ssum_LMS Ssum_LMS_norm ...
         Msum_LM Msum_LMS Msum_LMS_norm ...
         Lsum_LM Lsum_LMS Lsum_LMS_norm];
climW = [min(allW) max(allW)];

%% 2. Spatial maps: total surround weight from each cone class
% Columns: [LMS norm, LM, LMS]

figure;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
colormap(parula);

% ----- Row 1: S-cone surround weight -----
% (1,1) LMS_norm
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, Ssum_LMS_norm, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LMS norm: S surround total');
xlabel('deg'); ylabel('deg');

% (1,2) LM
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, Ssum_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LM: S surround total');
xlabel('deg'); ylabel('deg');

% (1,3) LMS
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, Ssum_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LMS: S surround total');
xlabel('deg'); ylabel('deg');

% ----- Row 2: M-cone surround weight -----
% (2,1) LMS_norm
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, Msum_LMS_norm, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LMS norm: M surround total');
xlabel('deg'); ylabel('deg');

% (2,2) LM
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, Msum_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LM: M surround total');
xlabel('deg'); ylabel('deg');

% (2,3) LMS
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, Msum_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LMS: M surround total');
xlabel('deg'); ylabel('deg');

% ----- Row 3: L-cone surround weight -----
% (3,1) LMS_norm
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, Lsum_LMS_norm, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LMS norm: L surround total');
xlabel('deg'); ylabel('deg');

% (3,2) LM
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, Lsum_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LM: L surround total');
xlabel('deg'); ylabel('deg');

% (3,3) LMS
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, Lsum_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim(climW);
title('LMS: L surround total');
xlabel('deg'); ylabel('deg');

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Total surround weight from that cone class';

%% 3. Cone-surround *counts* (|w| >= threshold), LM vs LMS vs LMS_norm
% Columns: [LMS norm, LM, LMS]

thr_LM  = mLM.minSurroundWeightForInclusionInComputing;
thr_LMS = mLMS.minSurroundWeightForInclusionInComputing;

% LM counts
numS_LM = sum(abs(W_S_LM) >= thr_LM, 1);
numM_LM = sum(abs(W_M_LM) >= thr_LM, 1);
numL_LM = sum(abs(W_L_LM) >= thr_LM, 1);

% LMS counts
numS_LMS = sum(abs(W_S_LMS) >= thr_LMS, 1);
numM_LMS = sum(abs(W_M_LMS) >= thr_LMS, 1);
numL_LMS = sum(abs(W_L_LMS) >= thr_LMS, 1);

% LMS_norm counts (use LMS threshold)
numS_LMS_norm = sum(abs(W_S_LMS_norm) >= thr_LMS, 1);
numM_LMS_norm = sum(abs(W_M_LMS_norm) >= thr_LMS, 1);
numL_LMS_norm = sum(abs(W_L_LMS_norm) >= thr_LMS, 1);

allCounts = [numS_LM numM_LM numL_LM ...
             numS_LMS numM_LMS numL_LMS ...
             numS_LMS_norm numM_LMS_norm numL_LMS_norm];
cminCnt = 0;
cmaxCnt = max(allCounts);
if cmaxCnt == 0, cmaxCnt = 1; end

figure;
tiledlayout(3,3,'TileSpacing','compact','Padding','compact');
colormap(parula);

% ----- Row 1: S-cone counts -----
% (1,1) LMS_norm
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numS_LMS_norm, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LMS norm: S surround count');
xlabel('deg'); ylabel('deg');

% (1,2) LM
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, numS_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LM: S surround count');
xlabel('deg'); ylabel('deg');

% (1,3) LMS
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numS_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LMS: S surround count');
xlabel('deg'); ylabel('deg');

% ----- Row 2: M-cone counts -----
% (2,1) LMS_norm
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numM_LMS_norm, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LMS norm: M surround count');
xlabel('deg'); ylabel('deg');

% (2,2) LM
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, numM_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LM: M surround count');
xlabel('deg'); ylabel('deg');

% (2,3) LMS
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numM_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LMS: M surround count');
xlabel('deg'); ylabel('deg');

% ----- Row 3: L-cone counts -----
% (3,1) LMS_norm
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numL_LMS_norm, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LMS norm: L surround count');
xlabel('deg'); ylabel('deg');

% (3,2) LM
nexttile;
scatter(pos_LM(:,1), pos_LM(:,2), 15, numL_LM, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LM: L surround count');
xlabel('deg'); ylabel('deg');

% (3,3) LMS
nexttile;
scatter(pos_LMS(:,1), pos_LMS(:,2), 15, numL_LMS, 'filled');
axis equal; xlim(xlims); ylim(ylims);
clim([cminCnt cmaxCnt]);
title('LMS: L surround count');
xlabel('deg'); ylabel('deg');

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = 'Number of cones in surround (w >= threshold)';

%% 4. Surround spectral tuning from wiring (LM vs LMS vs LMS_norm)

wave = cm_LM.wave(:);    % common wavelength sampling
qe   = cm_LM.qe;         % [nWaves x 3] (L,M,S)

qL = qe(:,1);
qM = qe(:,2);
qS = qe(:,3);

% LM model
specSurround_L_LM = qL * Lsum_LM;   % [nWaves x nRGC]
specSurround_M_LM = qM * Msum_LM;
specSurround_S_LM = qS * Ssum_LM;
specSurround_total_LM = specSurround_L_LM + ...
                        specSurround_M_LM + ...
                        specSurround_S_LM;

% LMS model
specSurround_L_LMS = qL * Lsum_LMS;
specSurround_M_LMS = qM * Msum_LMS;
specSurround_S_LMS = qS * Ssum_LMS;
specSurround_total_LMS = specSurround_L_LMS + ...
                         specSurround_M_LMS + ...
                         specSurround_S_LMS;

% LMS_norm model
specSurround_L_LMS_norm = qL * Lsum_LMS_norm;
specSurround_M_LMS_norm = qM * Msum_LMS_norm;
specSurround_S_LMS_norm = qS * Ssum_LMS_norm;
specSurround_total_LMS_norm = specSurround_L_LMS_norm + ...
                              specSurround_M_LMS_norm + ...
                              specSurround_S_LMS_norm;

%% 5. Overlaid surround spectral curves

meanSpecLM       = mean(specSurround_total_LM,       2);   % [nWaves x 1]
meanSpecLMS      = mean(specSurround_total_LMS,      2);
meanSpecLMS_norm = mean(specSurround_total_LMS_norm, 2);

yminSpec = min([specSurround_total_LM(:); ...
                specSurround_total_LMS(:); ...
                specSurround_total_LMS_norm(:)]);
ymaxSpec = max([specSurround_total_LM(:); ...
                specSurround_total_LMS(:); ...
                specSurround_total_LMS_norm(:)]);

figure; clf; hold on; set(gcf,'Color','w');

% All LM surrounds (light red)
plot(wave, specSurround_total_LM,       'Color', [1.0 0.7 0.7], 'LineWidth', 0.05);

% All LMS surrounds (light blue)
plot(wave, specSurround_total_LMS,      'Color', [0.7 0.7 1.0], 'LineWidth', 0.05);

% All LMS_norm surrounds (light cyan)
plot(wave, specSurround_total_LMS_norm, 'Color', [0.6 1.0 1.0], 'LineWidth', 0.05);

% Means: LM red, LMS blue, LMS_norm cyan
hLM       = plot(wave, meanSpecLM,       'r',                        'LineWidth', 2);
hLMS      = plot(wave, meanSpecLMS,      'b',                        'LineWidth', 2);
hLMS_norm = plot(wave, meanSpecLMS_norm, 'Color', [0.0 0.7 0.7],     'LineWidth', 2);

xlabel('Wavelength (nm)');
ylabel('Surround sensitivity (arb. units)');
title('mRGC surround spectral tuning: LM vs LMS vs LMS norm');
ylim([yminSpec ymaxSpec]);
legend([hLM hLMS hLMS_norm], {'LM mean','LMS mean','LMS norm mean'}, 'Location','best');
grid on;

%% 6. Paired dot plots: L, M, S surround totals (LM, LMS, LMS_norm) + differences

% Color convention (LM red, LMS blue, LMS_norm cyan)
LMcolor      = [0.8 0.1 0.1];   % LM
LMScolor     = [0.1 0.1 0.8];   % LMS
LMSnormColor = [0.0 0.6 0.6];   % LMS_norm
lineColor    = [0.8 0.8 0.8];   % connecting lines

% Jittered x positions so dots don't sit exactly on top of each other
baseNorm = 1;   % LMS_norm
baseLM   = 2;   % LM
baseLMS  = 3;   % LMS

offset = 0.04;
noise  = 0.04;

% M-cone surrounds
xM_norm = baseNorm - offset + noise*randn(1,nRGC);
xM_LM   = baseLM   + noise*randn(1,nRGC);
xM_LMS  = baseLMS  + offset + noise*randn(1,nRGC);

% L-cone surrounds
xL_norm = baseNorm - offset + noise*randn(1,nRGC);
xL_LM   = baseLM   + noise*randn(1,nRGC);
xL_LMS  = baseLMS  + offset + noise*randn(1,nRGC);

% S-cone surrounds
xS_norm = baseNorm - offset + noise*randn(1,nRGC);
xS_LM   = baseLM   + noise*randn(1,nRGC);
xS_LMS  = baseLMS  + offset + noise*randn(1,nRGC);

% Common y-limits from all L/M/S weights
allL = [Lsum_LM(:); Lsum_LMS(:); Lsum_LMS_norm(:)];
allM = [Msum_LM(:); Msum_LMS(:); Msum_LMS_norm(:)];
allS = [Ssum_LM(:); Ssum_LMS(:); Ssum_LMS_norm(:)];
ymin = min([allL; allM; allS]);
ymax = max([allL; allM; allS]);

% Per-RGC signed differences (LM vs LMS / LMS_norm)
Mdiff_orig = Msum_LM - Msum_LMS;        % LM - LMS
Ldiff_orig = Lsum_LM - Lsum_LMS;
Mdiff_norm = Msum_LM - Msum_LMS_norm;   % LM - LMS_norm
Ldiff_norm = Lsum_LM - Lsum_LMS_norm;

% LMS_norm - LMS (cone-class deltas)
Ldiff_normDelta = Lsum_LMS_norm - Lsum_LMS;   % LMS_norm - LMS (L)
Mdiff_normDelta = Msum_LMS_norm - Msum_LMS;   % LMS_norm - LMS (M)
Sdiff_normDelta = Ssum_LMS_norm - Ssum_LMS;   % LMS_norm - LMS (S)

% x-positions for LM-LMS vs LM-LMS_norm differences
baseM_diff      = 1;
baseM_norm_diff = 2;
baseL_diff      = 3;
baseL_norm_diff = 4;

diffOffset = 0.06;
diffNoise  = 0.03;

xDiff_M_orig = baseM_diff      - diffOffset + diffNoise*randn(1,nRGC);
xDiff_M_norm = baseM_norm_diff + diffOffset + diffNoise*randn(1,nRGC);

xDiff_L_orig = baseL_diff      - diffOffset + diffNoise*randn(1,nRGC);
xDiff_L_norm = baseL_norm_diff + diffOffset + diffNoise*randn(1,nRGC);

% x-positions for LMS_norm - LMS (L, M, S)
baseDiff_LMS_L = 1;
baseDiff_LMS_M = 2;
baseDiff_LMS_S = 3;

xDiff_L_normOnly = baseDiff_LMS_L - diffOffset + diffNoise*randn(1,nRGC);
xDiff_M_normOnly = baseDiff_LMS_M +            diffNoise*randn(1,nRGC);
xDiff_S_normOnly = baseDiff_LMS_S + diffOffset + diffNoise*randn(1,nRGC);

% Y limits across all difference measures
allDiffs = [Mdiff_orig(:); Ldiff_orig(:); ...
            Mdiff_norm(:); Ldiff_norm(:); ...
            Ldiff_normDelta(:); Mdiff_normDelta(:); Sdiff_normDelta(:)];
yDiffLim = [min(allDiffs) max(allDiffs)];

figure; clf;

%% TOP LEFT: M-cone surrounds (LMS norm, LM, LMS)
subplot(2,3,1); hold on;

for i = 1:nRGC
    % Lines between LMS_norm -> LM -> LMS
    plot([xM_norm(i) xM_LM(i)], [Msum_LMS_norm(i) Msum_LM(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
    plot([xM_LM(i)   xM_LMS(i)], [Msum_LM(i)      Msum_LMS(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
end

% LMS_norm (cyan rings)
plot(xM_norm, Msum_LMS_norm, 'o', ...
     'MarkerEdgeColor', LMSnormColor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

% LM (red rings)
plot(xM_LM,  Msum_LM, 'o', ...
     'MarkerEdgeColor', LMcolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

% LMS (blue rings)
plot(xM_LMS, Msum_LMS, 'o', ...
     'MarkerEdgeColor', LMScolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

xlim([0.5 3.5]);
ylim([ymin ymax]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'LMS norm','LM','LMS'});
ylabel('Total M-cone surround weight');
title('M surrounds');

%% TOP MIDDLE: L-cone surrounds (LMS norm, LM, LMS)
subplot(2,3,2); hold on;

for i = 1:nRGC
    plot([xL_norm(i) xL_LM(i)], [Lsum_LMS_norm(i) Lsum_LM(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
    plot([xL_LM(i)   xL_LMS(i)], [Lsum_LM(i)      Lsum_LMS(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
end

plot(xL_norm, Lsum_LMS_norm, 'o', ...
     'MarkerEdgeColor', LMSnormColor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

plot(xL_LM,  Lsum_LM, 'o', ...
     'MarkerEdgeColor', LMcolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

plot(xL_LMS, Lsum_LMS, 'o', ...
     'MarkerEdgeColor', LMScolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

xlim([0.5 3.5]);
ylim([ymin ymax]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'LMS norm','LM','LMS'});
ylabel('Total L-cone surround weight');
title('L surrounds');

%% TOP RIGHT: S-cone surrounds (LMS norm, LM, LMS)
subplot(2,3,3); hold on;

for i = 1:nRGC
    plot([xS_norm(i) xS_LM(i)], [Ssum_LMS_norm(i) Ssum_LM(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
    plot([xS_LM(i)   xS_LMS(i)], [Ssum_LM(i)      Ssum_LMS(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
end

plot(xS_norm, Ssum_LMS_norm, 'o', ...
     'MarkerEdgeColor', LMSnormColor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

plot(xS_LM,  Ssum_LM, 'o', ...
     'MarkerEdgeColor', LMcolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

plot(xS_LMS, Ssum_LMS, 'o', ...
     'MarkerEdgeColor', LMScolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

xlim([0.5 3.5]);
ylim([ymin ymax]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'LMS norm','LM','LMS'});
ylabel('Total S-cone surround weight');
title('S surrounds');

%% BOTTOM LEFT: LM - LMS differences (original and normalized)
subplot(2,3,4); hold on;
yline(0,'k--');

% Colors for differences
MdiffColor_orig = [0.0 0.5 0.0];    % darker green: LM - LMS (M)
MdiffColor_norm = [0.6 0.9 0.6];    % light green: LM - LMS_norm (M)

LdiffColor_orig = [0.4 0.0 0.4];    % darker purple: LM - LMS (L)
LdiffColor_norm = [0.85 0.7 0.85];  % light purple: LM - LMS_norm (L)

% M differences
hM_orig = plot(xDiff_M_orig, Mdiff_orig, 'o', ...
    'MarkerEdgeColor', MdiffColor_orig, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

hM_norm = plot(xDiff_M_norm, Mdiff_norm, 'o', ...
    'MarkerEdgeColor', MdiffColor_norm, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

% L differences
hL_orig = plot(xDiff_L_orig, Ldiff_orig, 'o', ...
    'MarkerEdgeColor', LdiffColor_orig, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

hL_norm = plot(xDiff_L_norm, Ldiff_norm, 'o', ...
    'MarkerEdgeColor', LdiffColor_norm, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

xlim([0.5 4.5]);
ylim(yDiffLim);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'M diff', 'norm M diff', 'L diff', 'norm L diff'});
ylabel('LM - LMS or LM - LMS norm surround weight');
title('Per-RGC LM vs LMS differences');

%% BOTTOM MIDDLE: LMS norm - LMS differences (L, M, S)
subplot(2,3,5); hold on;
yline(0,'k--');

SdiffColor_normOnly = [0.0 0.5 0.5];  % teal for S (norm - orig)

% L differences (norm - orig)
plot(xDiff_L_normOnly, Ldiff_normDelta, 'o', ...
     'MarkerEdgeColor', LdiffColor_norm, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

% M differences (norm - orig)
plot(xDiff_M_normOnly, Mdiff_normDelta, 'o', ...
     'MarkerEdgeColor', MdiffColor_norm, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

% S differences (norm - orig)
plot(xDiff_S_normOnly, Sdiff_normDelta, 'o', ...
     'MarkerEdgeColor', SdiffColor_normOnly, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

xlim([0.5 3.5]);
ylim(yDiffLim);
set(gca,'XTick',[1 2 3], 'XTickLabel', {'L diff','M diff','S diff'});
ylabel('LMS norm - LMS surround weight');
title('Per-RGC LMS norm - LMS');

%% LEFT: M-cone surrounds (LMS norm, LM, LMS)
subplot(1,3,1); hold on;

for i = 1:nRGC
    % Lines between LMS_norm -> LM -> LMS
    plot([xM_norm(i) xM_LM(i)], [Msum_LMS_norm(i) Msum_LM(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
    plot([xM_LM(i)   xM_LMS(i)], [Msum_LM(i)      Msum_LMS(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
end

% LMS_norm (cyan rings)
plot(xM_norm, Msum_LMS_norm, 'o', ...
     'MarkerEdgeColor', LMSnormColor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

% LM (red rings)
plot(xM_LM,  Msum_LM, 'o', ...
     'MarkerEdgeColor', LMcolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

% LMS (blue rings)
plot(xM_LMS, Msum_LMS, 'o', ...
     'MarkerEdgeColor', LMScolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

xlim([0.5 3.5]);
ylim([ymin ymax]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'LMS norm','LM','LMS'});
ylabel('Total M-cone surround weight');
title('M surrounds');

%% MIDDLE: L-cone surrounds (LMS norm, LM, LMS)
subplot(1,3,2); hold on;

for i = 1:nRGC
    plot([xL_norm(i) xL_LM(i)], [Lsum_LMS_norm(i) Lsum_LM(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
    plot([xL_LM(i)   xL_LMS(i)], [Lsum_LM(i)      Lsum_LMS(i)], '-', ...
         'Color', lineColor, 'LineWidth', 0.5);
end

plot(xL_norm, Lsum_LMS_norm, 'o', ...
     'MarkerEdgeColor', LMSnormColor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

plot(xL_LM,  Lsum_LM, 'o', ...
     'MarkerEdgeColor', LMcolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

plot(xL_LMS, Lsum_LMS, 'o', ...
     'MarkerEdgeColor', LMScolor, ...
     'MarkerFaceColor', 'none', ...
     'MarkerSize', 4);

xlim([0.5 3.5]);
ylim([ymin ymax]);
set(gca,'XTick',[1 2 3],'XTickLabel',{'LMS norm','LM','LMS'});
ylabel('Total L-cone surround weight');
title('L surrounds');

%% RIGHT: LM - LMS differences (original and normalized)
subplot(1,3,3); hold on;
yline(0,'k--');

% Colors for differences
MdiffColor_orig = [0.0 0.5 0.0];    % darker green: LM - LMS (M)
MdiffColor_norm = [0.6 0.9 0.6];    % light green: LM - LMS_norm (M)

LdiffColor_orig = [0.4 0.0 0.4];    % darker purple: LM - LMS (L)
LdiffColor_norm = [0.85 0.7 0.85];  % light purple: LM - LMS_norm (L)

% M differences
hM_orig = plot(xDiff_M_orig, Mdiff_orig, 'o', ...
    'MarkerEdgeColor', MdiffColor_orig, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

hM_norm = plot(xDiff_M_norm, Mdiff_norm, 'o', ...
    'MarkerEdgeColor', MdiffColor_norm, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

% L differences
hL_orig = plot(xDiff_L_orig, Ldiff_orig, 'o', ...
    'MarkerEdgeColor', LdiffColor_orig, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

hL_norm = plot(xDiff_L_norm, Ldiff_norm, 'o', ...
    'MarkerEdgeColor', LdiffColor_norm, ...
    'MarkerFaceColor', 'none', ...
    'MarkerSize', 4);

xlim([0.5 4.5]);
ylim(yDiffLim);
set(gca,'XTick',[1 2 3 4],'XTickLabel',{'M diff', 'norm M diff', 'L diff', 'norm L diff'});
ylabel('LM - LMS or norm LMS surround weight');
title('Per-RGC differences');

% legend([hM_orig hM_norm hL_orig hL_norm], ...
%        {'M: LM - LMS','M: LM - LMS norm', ...
%         'L: LM - LMS','L: LM - LMS norm'}, ...
%        'Location','best');
